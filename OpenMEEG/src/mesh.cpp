/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

The OpenMEEG software is a C++ package for solving the forward/inverse
problems of electroencephalography and magnetoencephalography.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's authors,  the holders of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#include <sstream>
#include <algorithm>
#include <mesh.h>
#include <Triangle_triangle_intersection.h>

namespace OpenMEEG {

    void Mesh::copy(const Mesh& m) {

        if (m.allocate_) {
            allocate_     = true;
            all_vertices_ = new Vertices; // allocates space for the vertices and copy
            all_vertices_->reserve(m.vertices().size()); 
            std::map<const Vertex*,Vertex*> map; // for the triangles
            for (const auto& vertex : m.vertices()) {
                add_vertex(*vertex);
                map[vertex] = vertices().back();
            }
            for (const auto& triangle : m.triangles()) {
                Triangle t(map[&triangle.vertex(0)],map[&triangle.vertex(1)],map[&triangle.vertex(2)]);
                triangles_.push_back(t);
            }
            update();
        } else {
            all_vertices_ = m.all_vertices_;
            set_vertices_ = m.set_vertices_;
            allocate_ = false;
            for (const auto& triangle : m.triangles())
                triangles_.push_back(triangle);
            build_mesh_vertices();
        }
        outermost_ = m.outermost_;
        name_      = m.name_;
    }

    /// Print informations about the mesh 

    void Mesh::info(const bool verbose) const {
        std::cout << "Info:: Mesh name/ID : "  << name() << std::endl;
        std::cout << "\t\t# vertices  : " << vertices().size() << std::endl;
        std::cout << "\t\t# triangles : " << triangles().size() << std::endl;
        std::cout << "\t\tEuler characteristic : " << vertices().size()-3*triangles().size()/2+triangles().size() << std::endl;

        double min_area = std::numeric_limits<double>::max();
        double max_area = 0.;
        for (const auto& triangle : triangles()) {
            min_area = std::min(triangle.area(),min_area);
            max_area = std::max(triangle.area(),max_area);
        }
        std::cout << "\t\tMin Area : " << min_area << std::endl;
        std::cout << "\t\tMax Area : " << max_area << std::endl;
        if (verbose) {
            std::cout << "Indices :" << std::endl;
            for (const auto& vertex : vertices())
                std::cout << "[" << *vertex << "] = " << vertex->index() << std::endl;
            for (const auto& triangle : triangles())
                std::cout << "[[" << triangle.vertex(0)
                          << "] , [" << triangle.vertex(1)
                          << "] , ["<< triangle.vertex(2)
                          << "]] \t = " << triangle.index() << std::endl;
        }
    }

    void Mesh::build_mesh_vertices() {

        // Sets do not preserve the order, and we would like to preserve it so we push_back in the vector as soon as the element is unique.

        std::set<const Vertex *> mesh_v;
        vertices_.clear();
        for (auto& triangle : triangles_)
            for (auto& vertex : triangle)
                if (mesh_v.insert(vertex).second)
                    vertices_.push_back(vertex);
    }

    void Mesh::destroy() {
        if (allocate_)
            delete all_vertices_;
        triangles_.clear();
        all_vertices_ = 0;
        vertices_.clear();
        set_vertices_.clear();
        name_.clear();
        links_.clear();
        outermost_ = false;
        allocate_ = false;
    }

    //  This is not the input of a mesh. It should be renamed...

    std::istream& operator>>(std::istream& is,Mesh& m) {
        unsigned vind[3];
        is >> vind[0] >> vind[1] >> vind[2];
        for (unsigned i=0;i<3;++i)
            if (vind[i]>=m.vertices_.size()) {
                std::cerr << "Unknown vertex: " << vind[i] << " (hint: vertex numbering often starts at 0). Aborting." << std::endl;
                exit(1);
            }
        Triangle t(m.vertices_[vind[0]],m.vertices_[vind[1]],m.vertices_[vind[2]]);
        m.triangles().push_back(t);
        return is;
    }

    /// properly add vertex to the list. (if not already added)
    void Mesh::add_vertex(const Vertex& v) {

        // try to insert the vertex to the set
        std::pair<std::set<Vertex>::iterator, bool> ret = set_vertices_.insert(v);
        if (ret.second) {
            // if inserted, then it is a new vertex, and we add it to both lists
            all_vertices_->push_back(v);
            vertices_.push_back(&(*all_vertices_->rbegin()));
        } else {
            // if not inserted, Either it belongs to another mesh or it was dupplicated in the same mesh
            // TODO this may take time for too big redundant mesh
            Vertices::iterator vit = std::find(all_vertices_->begin(), all_vertices_->end(), v);
            if (std::find(vertices_.begin(), vertices_.end(), &(*vit)) == vertices_.end()) {
                vertices_.push_back(&(*vit));
            }
        }
    }

    /// Update triangles area/normal, update links and vertices normals if needed
    void Mesh::update() {

        // empty unnecessary set

        set_vertices_.clear();

        make_adjacencies();

        // If indices are not set, we generate them for sorting edge and testing orientation
        if (allocate_) {
            if (vertices().front()->index()==unsigned(-1))
                generate_indices();
            correct_local_orientation();
        }

        // Compute triangles' normals and areas (after having the mesh locally reoriented)

        for (auto& triangle : triangles()) {
            Vect3 normaldir   = crossprod(triangle.vertex(0)-triangle.vertex(1),triangle.vertex(0)-triangle.vertex(2));
            triangle.area()   = normaldir.norm()/2.0;
            triangle.normal() = normaldir.normalize();
        }
    }

    /// Compute normals at vertices.

    Normal Mesh::normal(const Vertex& v) const {
        Normal _normal(0);
        for (VectPTriangle::const_iterator tit = links_.at(&v).begin(); tit != links_.at(&v).end(); ++tit)
            _normal += (*tit)->normal();
        _normal.normalize();
        return _normal;
    }

    /// Merge two meshes into one (without duplicating vertices).

    void Mesh::merge(const Mesh& m1, const Mesh& m2) {

        if (triangles_.size() != 0)
            warning("Mesh::merge Mesh must be empty.");

        allocate_ = true;
        all_vertices_ = new Vertices;
        all_vertices_->reserve(m1.vertices().size() + m2.vertices().size());

        std::map<unsigned,Vertex *> map1;
        for (const auto& vertex : m1.vertices()) {
            add_vertex(*vertex);
            map1[vertex->index()] = vertices().back();
        }

        std::map<unsigned,Vertex*> map2;
        for (const auto& vertex : m2.vertices()) {
            add_vertex(*vertex);
            map2[vertex->index()] = vertices().back();
        }

        for (const auto& triangle : m1.triangles())
            triangles_.push_back(Triangle(map1[triangle.vertex(0).index()], map1[triangle.vertex(1).index()], map1[triangle.vertex(2).index()]));

        for (const auto& triangle : m2.triangles())
            triangles_.push_back(Triangle(map2[triangle.vertex(0).index()], map2[triangle.vertex(1).index()], map2[triangle.vertex(2).index()]));

        update();
    }

    /// Smooth Mesh

    void Mesh::smooth(const double& smoothing_intensity,const unsigned& niter) {

        typedef std::vector<std::set<Vertex>> Neighbors;
        Neighbors neighbors(vertices().size());
        unsigned i = 0;
        for (const auto& vertex : vertices()) {
            for (const auto& triangle : adjacent_triangles(*vertex))
                for (unsigned  k=0;k<3;++k)
                    if (triangle->vertex(k)==*vertex) // TODO: Is a pointer comparison sufficient ?
                        neighbors[i].insert(triangle->vertex(k));
            ++i; // Can't we use the vertex index instead ? TODO
        }

        Vertices new_pts(vertices().size());
        for (unsigned n = 0; n < niter; ++n) {
            Neighbors::const_iterator nit = neighbors.begin();
            for (const auto& vertex : vertices()) {
                new_pts.push_back(*vertex);
                auto& newpt = new_pts.back();
                for (std::set<Vertex>::const_iterator it = nit->begin();it!=nit->end();++it)
                    newpt += (smoothing_intensity*(*it-*vertex))/nit->size();
                ++nit;
            }
            for (auto& vertex : vertices())
                *vertex = new_pts[i];
            new_pts.clear();
        }
        update(); // Updating triangles (areas + normals)
    }

    /// Sq. Norm Surface Gradient: square norm of the surfacic gradient of the P1 and P0 elements

    void Mesh::gradient_norm2(SymMatrix &A) const {

        /// Vertices (diagonal of the matrix).

        for (const auto& vp : vertices()) {
            const Vertex& v1 = *vp;
            const unsigned index = v1.index();
            for (const auto& t : links_.at(vp)) {
                const Edge& edge = t->edge(v1);
                const Vertex& v2 = edge.vertex(0);
                const Vertex& v3 = edge.vertex(1);
                A(index,index) += P1gradient(v1,v2,v3).norm2()*sqr(t->area());
            }
        }

        // edges

        for (const auto& triangle: triangles_) {
            const double area2 = sqr(triangle.area());
            const Edges& edges = triangle.edges();
            for (unsigned i=0;i<3;++i) {
                const Vertex& V1 = edges[i].vertex(0);
                const Vertex& V2 = edges[i].vertex(1);

                const unsigned ind1 = V1.index();
                const unsigned ind2 = V2.index();

                //  This is a symmetric matrix. We only need to fill the lower half matrix.

                if (ind1<ind2) {
                    const Vertex& V3 = triangle.vertex(i);
                    A(ind1,ind2) += dotprod(P1gradient(V1,V2,V3),P1gradient(V2,V3,V1))*area2;
                }
            }
        }

        // P0 gradients: loop on triangles
        if (!outermost_) // if it is an outermost mesh: p=0 thus no need for computing it
            for (const auto& triangle1 : triangles()) {
                A(triangle1.index(),triangle1.index()) = 0.;
                for (const auto& triangle2 : adjacent_triangles(triangle1))
                    if (triangle1.index()<triangle2->index()) // sym matrix only lower half
                        A(triangle1.index(), triangle2->index()) += P0gradient_norm2(triangle1,*triangle2)*triangle1.area()*triangle2->area();
            }
    }

    /// Laplacian Mesh: good approximation of Laplace-Beltrami operator
    // "Discrete Laplace Operator on Meshed Surfaces". by Belkin, Sun, Wang

    void Mesh::laplacian(SymMatrix& A) const {
        for (const auto& triangle : triangles())
            for (const auto& edge : triangle.edges()) {
                const Vertex& V1 = edge.vertex(0);
                const Vertex& V2 = edge.vertex(1);
                const unsigned ind1 = V1.index();
                const unsigned ind2 = V2.index();
                if (ind1<ind2) { // sym matrix only lower half
                    const double h = (V2-V1).norm();
                    A(ind1,ind2) += -triangle.area()/(12*Pi*sqr(h))*exp(-sqr(h)/(4*h));
                }
            }

        for (const auto& vertex : vertices())
            A(vertex->index(),vertex->index()) = -A.getlin(vertex->index()).sum();
    }

    bool Mesh::has_self_intersection() const {

        bool selfIntersects = false;
        for (auto tit1=triangles().begin();tit1!=triangles().end();++tit1)
            for (auto tit2=tit1;tit2!=triangles().end();++tit2)
                if (!tit1->contains(tit2->vertex(0)) && !tit1->contains(tit2->vertex(1)) && !tit1->contains(tit1->vertex(2)))
                    if (triangle_intersection(*tit1, *tit2)) {
                        selfIntersects = true;
                        std::cout << "Triangles " << tit1->index() << " and " << tit2->index() << " are intersecting." << std::endl;
                    }
        return selfIntersects;
    }

    double Mesh::solid_angle(const Vect3& p) const {
        double solangle = 0.0;
        for (const auto& triangle : triangles())
            solangle += p.solid_angle(triangle.vertex(0),triangle.vertex(1),triangle.vertex(2));
        return solangle;
    }

    bool Mesh::intersection(const Mesh& m) const {

        bool intersects = false;
        for (const auto& triangle1 : triangles())
            for (const auto& triangle2 : m.triangles())
                intersects = intersects | triangle_intersection(triangle1,triangle2);
        return intersects;
    }

    bool Mesh::triangle_intersection(const Triangle& T1, const Triangle& T2) const {

        const Vect3& p1 = T1.vertex(0);
        const Vect3& q1 = T1.vertex(1);
        const Vect3& r1 = T1.vertex(2);
        const Vect3& p2 = T2.vertex(0);
        const Vect3& q2 = T2.vertex(1);
        const Vect3& r2 = T2.vertex(2);

        double pp1[3] = {p1.x(), p1.y(), p1.z()};
        double qq1[3] = {q1.x(), q1.y(), q1.z()};
        double rr1[3] = {r1.x(), r1.y(), r1.z()};
        double pp2[3] = {p2.x(), p2.y(), p2.z()};
        double qq2[3] = {q2.x(), q2.y(), q2.z()};
        double rr2[3] = {r2.x(), r2.y(), r2.z()};
        return tri_tri_overlap_test_3d(pp1, qq1, rr1, pp2, qq2, rr2);
    }

    const Mesh::VectPTriangle& Mesh::get_triangles_for_vertex(const Vertex& V) const {

        std::map<const Vertex *, Mesh::VectPTriangle>::const_iterator it = links_.find(&V);
        if (it != links_.end()) {
            return it->second;
        } else {
            static Mesh::VectPTriangle a;
            return a;
        }
    }

    /// For IO:s -------------------------------------------------------------------------------------------

    unsigned Mesh::load(const std::string& filename, const bool& verbose, const bool& read_all) {

        if (triangles().size() != 0)
            destroy();

        if (read_all && (all_vertices_==0)) {
            unsigned nb_v = load(filename,false,false); // first allocates memory for the vertices
            all_vertices_ = new Vertices;
            all_vertices_->reserve(nb_v); 
            allocate_ = true;
        }

        std::string extension = getNameExtension(filename);
        std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
        unsigned return_value = 0;

        if (verbose)
            std::cout << "loading : " << filename << " as a \"" << extension << "\" file."<< std::endl;

        if (extension == std::string("vtk")) {
            return_value = load_vtk(filename, read_all);
        } else if (extension == std::string("tri")) {
            return_value = load_tri(filename, read_all);
        } else if (extension == std::string("bnd")) {
            return_value = load_bnd(filename, read_all);
        } else if (extension == std::string("mesh")) {
            return_value = load_mesh(filename, read_all);
        } else if (extension == std::string("off")) {
            return_value = load_off(filename, read_all);
        } else if (extension == std::string("gii")) {
            return_value = load_gifti(filename, read_all);
        } else {
            std::cerr << "IO: load: Unknown mesh file format for " << filename << std::endl;
            exit(1);
        }

        if (read_all)
            update();

        if (verbose)
            info();

        if (allocate_ && read_all) // we generates the indices of these mesh vertices
            generate_indices();

        return return_value;
    }

    void Mesh::generate_indices() {
        unsigned index = 0;
        for (auto& vertex : vertices())
            vertex->index() = index++;
        for (auto& triangle : triangles())
            triangle.index() = index++;
    }

    void Mesh::save(const std::string& filename) const {

        std::string extension = getNameExtension(filename);

        std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);

        if (extension==std::string("vtk")) {
            save_vtk(filename);
        } else if (extension==std::string("tri")) {
            save_tri(filename);
        } else if (extension==std::string("bnd")) {
            save_bnd(filename);
        } else if (extension==std::string("mesh")) {
            save_mesh(filename);
        } else if (extension==std::string("off")) {
            save_off(filename);
        } else if (extension==std::string("gii")) {
            save_gifti(filename);
        } else {
            std::cerr << "Unknown file format for : " << filename << std::endl;
            exit(1);
        }
    }

    #ifdef USE_VTK
    unsigned Mesh::get_data_from_vtk_reader(vtkPolyDataReader* reader, const bool& read_all) {

        reader->Update();
        vtkPolyData *vtkMesh = reader->GetOutput();

        unsigned npts, ntrgs;
        npts = vtkMesh->GetNumberOfPoints();

        if (!read_all)
            return npts;

        for (unsigned i=0;i<npts;++i)
            add_vertex(Vertex(vtkMesh->GetPoint(i)[0], vtkMesh->GetPoint(i)[1], vtkMesh->GetPoint(i)[2]));

        ntrgs = vtkMesh->GetNumberOfCells();

        vtkIdList *l;
        triangles().reserve(ntrgs);

        for (unsigned i = 0; i < ntrgs; ++i)
            if (vtkMesh->GetCellType(i) == VTK_TRIANGLE) {
                l = vtkMesh->GetCell(i)->GetPointIds();
                triangles().push_back(Triangle(vertices()[l->GetId(0)],vertices()[l->GetId(1)],vertices()[l->GetId(2)]));
            } else {
                std::cerr << "Mesh \"" << name_ << "\" is not a triangulation" << std::endl;
                exit(1);
            }

        return 0;
    }

    unsigned Mesh::load_vtk(std::istream& is, const bool& read_all) {

        // get length of file:
        is.seekg (0, ios::end);
        int length = is.tellg();
        is.seekg (0, ios::beg);

        // allocate memory:
        char * buffer = new char [length];

        // read data as a block:
        is.read (buffer, length);

        // held buffer by the array buf:
        vtkCharArray* buf = vtkCharArray::New();
        buf->SetArray(buffer, length, 1);

        vtkPolyDataReader* reader = vtkPolyDataReader::New();
        reader->SetInputArray(buf); // Specify 'buf' to be used when reading from a string
        reader->SetReadFromInputString(1);  // Enable reading from the InputArray 'buf' instead of the default, a file

        unsigned return_value = 0;
        return_value = get_data_from_vtk_reader(reader, read_all);

        delete[] buffer;
        reader->Delete();

        return return_value;
    }

    unsigned Mesh::load_vtk(const std::string& filename, const bool& read_all) {

        std::string s = filename;
        vtkPolyDataReader *reader = vtkPolyDataReader::New();
        reader->SetFileName(filename.c_str()); // Specify file name of vtk data file to read
        if (!reader->IsFilePolyData()) {
            std::cerr << "Mesh \"" << name_ << "\" is not a valid vtk poly data file" << std::endl;
            reader->Delete();
            exit(1);
        }

        unsigned return_value = 0;
        return_value = get_data_from_vtk_reader(reader, read_all);
        return return_value;
    }
    #endif

    #ifdef USE_GIFTI
    unsigned Mesh::load_gifti(const std::string& filename, const bool& read_all) {   

        // Use gifti_io API
        int read_data = 0; 
        gifti_image* gim = gifti_read_image(filename.c_str(), read_data);

        if (gim->numDA!=2)
            throw std::invalid_argument("OpenMEEG only handles gifti files containing two arrays.");

        // find which array contains the points and which the triangles
        unsigned ipts, itrgs;
        unsigned iit = 0;
        while (iit < gim->numDA) {
            if (gim->darray[iit]->intent == NIFTI_INTENT_POINTSET) {
                ipts = iit;
            } else if (gim->darray[iit]->intent == NIFTI_INTENT_TRIANGLE) {
                itrgs = iit;
            }
            ++iit;
        }
        if ((gim->darray[ipts]->dims[1]!=3) || // 3D points
            (gim->darray[itrgs]->dims[1]!=3))  // 3 indices per triangle
            throw std::invalid_argument("OpenMEEG only handles 3D surfacic meshes.");

        unsigned npts  = gim->darray[ipts]->dims[0];
        unsigned ntrgs = gim->darray[itrgs]->dims[0];

        gifti_free_image(gim);
        
        if (!read_all) { 
            return npts; 
        }
        read_data = 1; // now, load the data
        gim = gifti_read_image(filename.c_str(), read_data);
        
        float * pts_data;
        if (gim->darray[ipts]->datatype == NIFTI_TYPE_FLOAT32 ) { // float
            pts_data = (float *)gim->darray[ipts]->data;
        // } else if (gim->darray[ipts]->datatype == NIFTI_TYPE_FLOAT64 ) { // double
            // pts_data = (double *)gim->darray[ipts]->data;
        } else {
            throw std::invalid_argument("OpenMEEG:: gifti mesh points not float nor double !?.");
        }
        if (gim->darray[ipts]->ind_ord == GIFTI_IND_ORD_ROW_MAJOR ) { // RowMajor
            for (unsigned i = 0; i < npts; ++i) {
                add_vertex(Vertex(pts_data[i*3],pts_data[i*3+1],pts_data[i*3+2]));
            }
        } else {
            for (unsigned i = 0; i < npts; ++i) {
                add_vertex(Vertex(pts_data[i],pts_data[i+npts],pts_data[i+2*npts]));
            }
        }

        reserve(ntrgs);
        int * trgs = (int *)gim->darray[itrgs]->data;

        if (gim->darray[itrgs]->ind_ord == GIFTI_IND_ORD_ROW_MAJOR ) { // RowMajor
            for (unsigned i = 0; i < ntrgs; ++i) {
                push_back(Triangle(vertices_[trgs[i*3]],vertices_[trgs[i*3+1]],vertices_[trgs[i*3+2]]));
            }
        } else {
            for (unsigned i = 0; i < ntrgs; ++i) {
                push_back(Triangle(vertices_[trgs[i]],vertices_[trgs[i+ntrgs]],vertices_[trgs[i+2*ntrgs]]));
            }
        }

        // free all
        gifti_free_image(gim);
        return 0;
    }
    #endif

    unsigned Mesh::load_mesh(std::istream& is, const bool& read_all) {

        unsigned char* uc = new unsigned char[5]; // File format
        is.read((char*)uc, sizeof(unsigned char)*5);
        delete[] uc;

        uc = new unsigned char[4]; // lbindian
        is.read((char*)uc, sizeof(unsigned char)*4);
        delete[] uc;

        unsigned int* ui = new unsigned int[1]; // arg_size
        is.read((char*)ui, sizeof(unsigned int));
        unsigned int arg_size = ui[0];
        delete[] ui;

        uc = new unsigned char[arg_size]; // Trash
        is.read((char*)uc, sizeof(unsigned char)*arg_size);
        delete[] uc;

        ui = new unsigned int[1]; // vertex_per_face
        is.read((char*)ui, sizeof(unsigned int));
        unsigned int vertex_per_face = ui[0];
        delete[] ui;

        ui = new unsigned int[1]; // mesh_time
        is.read((char*)ui, sizeof(unsigned int));
        unsigned int mesh_time = ui[0];
        delete[] ui;

        ui = new unsigned int[1]; // mesh_step
        is.read((char*)ui, sizeof(unsigned int));
        delete[] ui;

        ui = new unsigned int[1]; // vertex number
        is.read((char*)ui, sizeof(unsigned int));
        unsigned npts = ui[0];
        delete[] ui;
        
        if (!read_all)
            return npts; 

        if (vertex_per_face!=3) // Support only for triangulations
            throw std::invalid_argument("OpenMEEG only handles 3D surfacic meshes.");
        if (mesh_time!=1) // Support only 1 time frame
            throw std::invalid_argument("OpenMEEG only handles 3D surfacic meshes with one time frame.");

        float* pts_raw = new float[npts*3]; // Points
        is.read((char*)pts_raw, sizeof(float)*npts*3);

        ui = new unsigned int[1]; // arg_size
        is.read((char*)ui, sizeof(unsigned int));
        delete[] ui;

        float* normals_raw = new float[npts*3]; // Normals
        is.read((char*)normals_raw, sizeof(float)*npts*3);

        ui = new unsigned int[1]; // arg_size
        is.read((char*)ui, sizeof(unsigned int));
        delete[] ui;

        ui = new unsigned int[1]; // number of faces
        is.read((char*)ui, sizeof(unsigned int));
        unsigned ntrgs;
        ntrgs = ui[0];
        delete[] ui;

        unsigned int* faces_raw = new unsigned int[ntrgs*3]; // Faces
        is.read((char*)faces_raw, sizeof(unsigned int)*ntrgs*3);

        for (unsigned i = 0; i < npts; ++i)
            add_vertex(Vertex(pts_raw[i*3+0],pts_raw[i*3+1],pts_raw[i*3+2]));

        triangles().reserve(ntrgs);

        for (unsigned i=0;i<ntrgs;++i)
            triangles().push_back(Triangle(vertices_[faces_raw[i*3+0]],vertices_[faces_raw[i*3+1]],vertices_[faces_raw[i*3+2]]));

        delete[] faces_raw;
        delete[] normals_raw;
        delete[] pts_raw;
        return 0;
    }

    unsigned Mesh::load_mesh(const std::string& filename, const bool& read_all) {

        std::ifstream f(filename.c_str(), std::ios::binary);
        if (!f.is_open()) {
            std::ostringstream ost;
            ost << "Error opening MESH file: " << filename << std::endl;
            throw std::invalid_argument(ost.str());
        }
        unsigned return_value = 0;
        return_value = load_mesh(f, read_all);
        f.close();
        return return_value;
    }

    unsigned Mesh::load_tri(std::istream& f,const bool& read_all) {

        f.seekg(0,std::ios_base::beg);

        char ch;
        unsigned npts, ntrgs;
        f >> ch;
        f >> npts;

        if (!read_all)
            return npts;

        for (unsigned i = 0; i < npts; ++i) {
            Vertex v;
            Normal n;
            f >> v >> n;
            add_vertex(v);
        }
        f >> ch >> ntrgs >> ntrgs >> ntrgs; // This number is repeated 3 times

        triangles().reserve(ntrgs);
        for (unsigned i=0;i<ntrgs;++i)
            f >> *this;

        return ntrgs;
    }

    unsigned Mesh::load_tri(const std::string& filename, const bool& read_all) {

        std::string s = filename;
        std::ifstream f(filename.c_str());
        if (!f.is_open()) {
            std::ostringstream ost;
            ost << "Error opening TRI file: " << filename << std::endl;
            throw std::invalid_argument(ost.str());
        }
        unsigned return_value = 0;
        return_value = load_tri(f, read_all);
        f.close();

        return return_value;
    }

    unsigned Mesh::load_bnd(std::istream& f, const bool& read_all) {

        std::string line;
        std::string st;

        f.seekg( 0, std::ios_base::beg);

        f >> io_utils::skip_comments('#') >> st;
        if (st == "Type=") {
            io_utils::skip_line(f);
            f >> io_utils::skip_comments('#') >> st;
        }

        om_error(st == "NumberPositions=");
        unsigned npts, ntrgs;
        f >> npts;

        if (!read_all)
            return npts;

        f >> io_utils::skip_comments('#') >> st;
        if (st == "UnitPosition")
            io_utils::skip_line(f); // skip : "UnitPosition mm"

        f >> io_utils::skip_comments('#') >> st;
        om_error(st == "Positions");

        for( unsigned i = 0; i < npts; ++i) {
            Vertex v;
            f >> io_utils::skip_comments('#') >> v;
            add_vertex(v);
        }

        f >> io_utils::skip_comments('#') >> st;
        om_error(st == "NumberPolygons=");
        f >> io_utils::skip_comments('#') >> ntrgs;

        f >> io_utils::skip_comments('#') >> st;
        om_error(st == "TypePolygons=");
        f >> io_utils::skip_comments('#') >> st;
        om_error(st == "3");

        f >> io_utils::skip_comments('#') >> st;
        om_error(st == "Polygons");

        triangles().reserve(ntrgs);

        for (unsigned i=0; i<ntrgs; ++i)
             f >> io_utils::skip_comments('#') >> *this;

        return 0;
    }

    unsigned Mesh::load_bnd(const std::string& filename, const bool& read_all) {        

        std::string s = filename;
        std::ifstream f(filename.c_str());

        if (!f.is_open()) {
            std::cerr << "Error opening BND file: " << filename << std::endl;
            exit(1);
        }
        unsigned return_value = 0;
        return_value = load_bnd(f, read_all);
        f.close();

        return return_value;
    }

    unsigned Mesh::load_off(std::istream& f, const bool& read_all) {
        char tmp[128];
        int trash;
        f >> tmp;        // put the "OFF" string
        unsigned npts, ntrgs;
        f >> npts;
        f >> ntrgs;
        f >> trash;

        if (!read_all)
            return npts;

        for (unsigned i=0; i<npts; ++i) {
            Vertex v;
            f >> v;
            add_vertex(v);
        }

        triangles().reserve(ntrgs);

        for (unsigned i=0; i<ntrgs; ++i)
            f >> trash >> *this;        // put the "3" to trash

        return 0;
    }

    unsigned Mesh::load_off(const std::string& filename, const bool& read_all) {

        std::string s = filename;
        std::ifstream f(filename.c_str());
        if (!f.is_open()) {
            std::cerr << "Error opening OFF file: " << filename << std::endl;
            exit(1);
        }
        unsigned return_value = 0;
        return_value = load_off(f, read_all);
        f.close();
        return return_value;
    }
    
    void Mesh::save_vtk(const std::string& filename) const {

        std::ofstream os(filename.c_str());
        os << "# vtk DataFile Version 2.0" << std::endl;
        os << "File " << filename << " generated by OpenMEEG" << std::endl;
        os << "ASCII" << std::endl;
        os << "DATASET POLYDATA" << std::endl;
        os << "POINTS " << vertices().size() << " float" << std::endl;

        std::map<const Vertex *, unsigned> map;
        unsigned i = 0;
        for (const auto& vertex : vertices()) {
            map[vertex] = i++;
            os << *vertex << std::endl;
        }
        os << "POLYGONS " << triangles().size() << " " << triangles().size()*4 << std::endl;
        for (const auto& triangle : triangles())
            os << "3 " << map[&(triangle.vertex(0))] << " " << map[&(triangle.vertex(1))] << " " << map[&(triangle.vertex(2))] << std::endl;

        os << "CELL_DATA " << triangles().size() << std::endl;
        os << "POINT_DATA " << vertices().size() << std::endl;
        os << "NORMALS normals float" << std::endl;
        for (const auto& vertex : vertices())
            os << normal(*vertex) << std::endl;

        os.close();
    }

    void Mesh::save_bnd(const std::string& filename) const {

        std::ofstream os(filename.c_str());
        os << "# Bnd mesh file generated by OpenMeeg" << std::endl;
        os << "Type= Unknown" << std::endl;
        os << "NumberPositions= " << vertices().size() << std::endl;
        os << "UnitPosition\tmm" << std::endl;
        os << "Positions" << std::endl;
        std::map<const Vertex *, unsigned> map;
        unsigned i = 0;
        for (const auto& vertex : vertices()) {
            map[vertex] = i++;
            os << *vertex << std::endl;
        }
        os << "NumberPolygons= " << triangles().size() << std::endl;
        os << "TypePolygons=\t3" << std::endl;
        os << "Polygons" << std::endl;
        for (const auto& triangle : triangles())
            os << map[&(triangle.vertex(0))] << " " << map[&(triangle.vertex(1))] << " " << map[&(triangle.vertex(2))] << std::endl;

        os.close();
    }

    void Mesh::save_tri(const std::string& filename) const {

        std::ofstream os(filename.c_str());
        os << "- " << vertices().size() << std::endl;
        std::map<const Vertex *, unsigned> map;
        unsigned i = 0;
        for (const auto& vertex : vertices()) {
            map[vertex] = i++;
            os << *vertex << " " << normal(*vertex) << std::endl;
        }
        os << "- " << triangles().size() << " " << triangles().size() << " " << triangles().size() << std::endl;
        for (const auto& triangle : triangles())
            os << map[&(triangle.vertex(0))] << " " << map[&(triangle.vertex(1))] << " " << map[&(triangle.vertex(2))] << std::endl;

        os.close();
    }

    void Mesh::save_off(const std::string& filename) const {

        std::ofstream os(filename.c_str());
        os << "OFF" << std::endl;
        os << vertices().size() << " " << triangles().size() << " 0" << std::endl;
        std::map<const Vertex *, unsigned> map;

        unsigned i = 0;
        for (const auto& vertex : vertices()) {
            map[vertex] = i++;
            os << *vertex << std::endl;
        }

        for (const auto& triangle : triangles())
            os << "3 " << map[&(triangle.vertex(0))] << " " << map[&(triangle.vertex(1))] << " " << map[&(triangle.vertex(2))] << std::endl;

        os.close();
    }

    void Mesh::save_mesh(const std::string& filename) const {

        std::ofstream os(filename.c_str(), std::ios::binary);

        unsigned char format[5] = {'b', 'i', 'n', 'a', 'r'}; // File format
        os.write((char*)format, sizeof(unsigned char)*5);

        unsigned char lbindian[4] = {'D', 'C', 'B', 'A'}; // lbindian
        os.write((char*)lbindian, sizeof(unsigned char)*4);

        unsigned int arg_size[1] = {4}; // arg_size
        os.write((char*)arg_size, sizeof(unsigned int));

        unsigned char VOID[4] = {'V', 'O', 'I', 'D'}; // Trash
        os.write((char*)VOID, sizeof(unsigned char)*4);

        unsigned int vertex_per_face[1] = {3}; // vertex_per_face
        os.write((char*)vertex_per_face, sizeof(unsigned int));

        unsigned int mesh_time[1] = {1}; // mesh_time
        os.write((char*)mesh_time, sizeof(unsigned int));

        unsigned int mesh_step[1] = {0}; // mesh_step
        os.write((char*)mesh_step, sizeof(unsigned int));

        unsigned int vertex_number[1] = {static_cast<unsigned>(vertices().size())}; // vertex number
        os.write((char*)vertex_number, sizeof(unsigned int));

        float* pts_raw = new float[vertices().size()*3]; // Points
        float* normals_raw = new float[vertices().size()*3]; // Normals
        unsigned int* faces_raw = new unsigned int[triangles().size()*3]; // Faces

        std::map<const Vertex *, unsigned> map;
        unsigned i = 0;
        for (const auto& vertex : vertices()) {
            map[vertex] = i;
            pts_raw[i*3+0]     = static_cast<float>(vertex->x());
            pts_raw[i*3+1]     = static_cast<float>(vertex->y());
            pts_raw[i*3+2]     = static_cast<float>(vertex->z());
            const Normal& n = normal(*vertex);
            normals_raw[i*3+0] = static_cast<float>(n.x());
            normals_raw[i*3+1] = static_cast<float>(n.y());
            normals_raw[i*3+2] = static_cast<float>(n.z());
            ++i;
        }

        i = 0;
        for (const auto& triangle : triangles()) {
            faces_raw[i*3+0] = map[&(triangle.vertex(0))];
            faces_raw[i*3+1] = map[&(triangle.vertex(1))];
            faces_raw[i*3+2] = map[&(triangle.vertex(2))];
            ++i;
        }

        os.write((char*)pts_raw, sizeof(float)*vertices().size()*3);           // vertices
        os.write((char*)vertex_number, sizeof(unsigned int));              // arg size : npts
        os.write((char*)normals_raw, sizeof(float)*vertices().size()*3);       // normals
        unsigned char zero[1] = {0};
        os.write((char*)zero, sizeof(char));                               // arg size : 0
        unsigned int faces_number[1] = {static_cast<unsigned>(triangles().size())};
        os.write((char*)faces_number, sizeof(unsigned int));               // ntrgs
        os.write((char*)faces_raw, sizeof(unsigned int)*triangles().size()*3); // triangles

        delete[] faces_raw;
        delete[] normals_raw;
        delete[] pts_raw;
        os.close();
    }

    #ifdef USE_GIFTI
    void Mesh::save_gifti(const std::string& filename) const {

        int dims[2] = {int(vertices().size()), 3};
        gifti_image * gim = gifti_create_image(1, NIFTI_INTENT_POINTSET, NIFTI_TYPE_FLOAT32, 2, dims, 1);
        gim->darray[0]->encoding = GIFTI_ENCODING_B64GZ;

        // vertices
        std::map<const Vertex *, int> mapVertexIndex;
        gim->darray[0]->data = new float[vertices().size()*3];
        float * pts = (float *)gim->darray[0]->data;
        int i = 0;
        for (const_vertex_iterator vit = vertex_begin(); vit != vertex_end(); ++vit, ++i) {
            mapVertexIndex[*vit] = i;
            pts[3*i]   = (*vit)->x();
            pts[3*i+1] = (*vit)->y();
            pts[3*i+2] = (*vit)->z();
        }

        // triangles
        gifti_add_empty_darray(gim, 1);
        giiDataArray *gDA = gim->darray[1];
        gifti_set_DA_defaults(gDA);
        gDA->intent   = NIFTI_INTENT_TRIANGLE;
        gDA->datatype = NIFTI_TYPE_INT32;
        gDA->ind_ord  = GIFTI_IND_ORD_ROW_MAJOR;
        gDA->num_dim  = 2;
        gDA->dims[0]  = int(triangles().size());
        gDA->dims[1]  = 3;
        gDA->nvals    = 3*triangles().size();
        gDA->encoding = GIFTI_ENCODING_B64GZ;
        gDA->endian   = GIFTI_ENDIAN_LITTLE;

        gDA->data = new int[triangles().size()*3];
        int * trgs = (int *)gDA->data;
        i = 0;
        for (const_iterator tit = begin(); tit != end(); ++tit, ++i) {
            trgs[3*i]   = mapVertexIndex[&(tit->vertex(0))];
            trgs[3*i+1] = mapVertexIndex[&(tit->vertex(1))];
            trgs[3*i+2] = mapVertexIndex[&(tit->vertex(2))];
        }
        // gifti_add_to_meta(&gDA->meta, "TopologicalType", "Closed", 0);
        gifti_add_to_meta(&gDA->meta, "Name",name_.c_str(), 0);

        gifti_write_image(gim, filename.c_str(), 1);
    }
    #endif

    const Mesh::EdgeMap Mesh::compute_edge_map() const {

        // Associate an integer with each edge.
        // Well oriented inner edges will be mapped to 0, border edges to 1 and badly oriented edges to 2.
        // The algorithm goes through each triangle edge e=(first vertex, second vertex)
        // If e is ordered with (lower index, higher index) add 1 to its map else remove 1.

        EdgeMap edgemap;
        auto lambda = [&](const Vertex& V1,const Vertex& V2,const int incr) {
            const auto& index = std::make_pair(&V1,&V2);
            if (edgemap.count(index)==0)
                edgemap[index] = incr;
            else
                edgemap[index] += incr;
        };

        for (const auto& triangle : triangles())
            for (const auto& edge : triangle.edges()) {
                const Vertex& V1 = edge.vertex(0);
                const Vertex& V2 = edge.vertex(1);
                if (V1.index()>V2.index())
                    lambda(V1,V2,1);
                else
                    lambda(V2,V1,-1);
            }

        return edgemap;
    }

    void Mesh::correct_local_orientation() {
        if (!has_correct_orientation()) {
            std::cerr << "Reorienting..." << std::endl << std::endl;
            std::stack<Triangle*>    triangle_stack;
            std::map<Triangle*,bool> reoriented_triangles;
            triangle_stack.push(&triangles().front());
            reoriented_triangles[&triangles().front()] = true;

            const auto has_same_edge = [](const Edges& edges1,const Edges& edges2) {
                for (const auto& edge2 : edges2)
                    for (const auto& edge1 : edges1)
                        if (edge1==edge2)
                            return true;
                return false;
            };

            while (!triangle_stack.empty()) {
                const Triangle& t1     = *triangle_stack.top();
                const Edges&    edges1 = t1.edges();
                triangle_stack.pop();
                for (const auto& tp : adjacent_triangles(t1))
                    if (reoriented_triangles.count(tp)==0) {
                        triangle_stack.push(tp);
                        Triangle& t2 = *tp;
                        const Edges& edges2 = t2.edges();
                        if (has_same_edge(edges1,edges2))
                            t2.change_orientation();
                        reoriented_triangles[tp] = true;
                    }
            }
        }
    }

    /// warning: a mesh may not be closed (as opposite to an interface)
    /// this function is mainly needed when meshing closed mesh with CGAL.

    void Mesh::correct_global_orientation() {

        Vect3 center(0.,0.,0.);
        for (const auto& vertex : vertices())
            center += *vertex;

        center /= vertices().size();
        const double solangle = solid_angle(center);
        if (!almost_equal(solangle,-4*Pi)) {
            if (almost_equal(solangle,0.0)) {
                std::cout << "Center point :" << center << " is on the mesh." << std::endl;
            } else if (almost_equal(solangle,4*Pi)) {
                change_orientation();
            } else {
                std::cout << "Not a closed mesh." << std::endl;
            }
        }
    }

    bool Mesh::has_correct_orientation() const {

        /// Check the local orientation (that all the triangles are all oriented in the same way)

        const EdgeMap mape = compute_edge_map();

        for (EdgeMap::const_iterator eit = mape.begin(); eit != mape.end(); ++eit)
            if (std::abs(eit->second) == 2) {
                std::cerr << "Local orientation problem..." << std::endl << std::endl;
                return false;
            }

        return true;
    }
}
