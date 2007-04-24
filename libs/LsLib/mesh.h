/*! \file
    \brief file containing the mesh class.
*/
#ifndef MESH_H
#define MESH_H

#include <set>
#include "point.h"

#if WIN32
#   pragma warning (disable: 4702)
#   pragma warning (disable: 4127) // for : mesh.h(367) : warning C4127: conditional expression is constant
#endif

//! set of int described using STL.
typedef std::set<int> intSet;

#ifdef USE_VTK
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyDataReader.h>
#include <vtkCellArray.h>
#include <vtkProperty.h>
#include <vtkPolyDataNormals.h>
#include <vtkCleanPolyData.h>
#endif

//! Different sort of mesh file format.
enum MeshFileType {NotSpecified,MFT_VTK,MFT_3D,MFT_OFF,MFT_MESH,MFT_TRI,MFT_MD,MFT_DSGL};

/*! Template for mesh.
    \param dim dimension of the points composing the mesh (2 for 2D, 3 for 3D, etc.).
    \attention mesh (lslib) is different of Mesh (winlib).
*/
template <int dim> class mesh
{
    /*! \var npts Number of points in the mesh.
    */
    int npts;
    
    /*! \var ncls Number of cells in the mesh.
    */
    int ncls;
    
    rPoint(dim) *pts; //!< Points composing the mesh.
    bool withNormals; //!< Indicate the presence of normals.
    /*! \brief Normals of the mesh used for display.
        \attention There is only one normal for each point and not one for each point for each cell.
    */
    rPoint(dim) *normals;

    iPoint(dim) *cells; //!< Cells composing the mesh.
    int* count; //!< Counter used to allow that the copies of same mesh share the same data.

    /*! Subfunction of the destructor of mesh.
        \attention Destroy shared data only if it is the last copy that uses those data.
    */
    void kill()
    {
        if (!count)
            return;
        (*count)--;
        if (!(*count)) {
            delete[] pts;
            delete[] cells;
            if (withNormals)
                delete[] normals;
            delete count;
        }
    }

    /*! Allocate the memory for the mesh.
        \param np Number of points.
        \param nc Number of cells.
        \param wn
            - \c true with normals.
            - \c false without normals.
    */
    void alloc(int np,int nc,bool wn) {
        npts=np;
        pts=new rPoint(dim)[np];
        ncls=nc;
        cells=new iPoint(dim)[nc];
        withNormals=wn;
        if (withNormals)
            normals=new rPoint(dim)[np];
        count=new int(1);
    }

    /*! Copy a mesh.
        \param M Mesh to copy.
        \attention No new memory is allocated, this copy is a link to the source mesh.
    */
    void copy(const mesh<dim>&M)
    {
        if (!M.count) {
            count=0;
            return;
        }
        npts=M.npts;
        ncls=M.ncls;
        pts=M.pts;
        cells=M.cells;
        count=M.count;
        withNormals=M.withNormals;
        if (withNormals)
            normals=M.normals;
        (*count)++;
    }

public:
    /*! Clone a mesh.
        \return a copy of the mesh.
        \attention This is not a link to the mesh. New memory space is allocated and data are copied.
    */
    mesh<dim> clone() const
    {
        if (!count)
            return mesh();
        mesh<dim> M(npts,ncls,withNormals);
        for (int i=0;i<npts;i++)
            M.pts[i]=pts[i];
        if (withNormals)
            for (int i=0;i<npts;i++)
                M.normals[i]=normals[i];
        for (int i=0;i<ncls;i++)
            M.cells[i]=cells[i];
        return M;
    }

    //! Constructor.
    mesh(void)
    {
        count=0;
    }
    //! Destructor.
    ~mesh(void)
    {
        kill();
    }
    /*! Constructor and initializer.
        \param np Number of points.
        \param nc Number of cells.
        \param wn "with normals"
    */
    mesh(int np,int nc,bool wn)
    {
        alloc(np,nc,wn);
    }
    //! Copy constructor.
    mesh(const mesh<dim>& M)
    {
        copy(M);
    }
    /*! Create a new mesh with different dimension.
        \param np Number of points.
        \param nc Number of cells.
        \attention Data are lost.
    */
    void resize(int np,int nc)
    {
        kill();
        alloc(np,nc,withNormals);
    }

    //! Operator= .
    const mesh<dim>& operator=(const mesh<dim>& M)
    {
        kill();
        copy(M);
        return *this;
    }

    inline int npoints() const { return npts; } //!< Accessor to the number of points.
    inline int ncells() const { return ncls; } //!< Accessor to the number of cells.
    inline rPoint(dim) point(int i) const { return pts[i]; } //!< Accessor to a point (read).
    inline rPoint(dim)& point(int i) { return pts[i]; } //!< Accessor to a point (read/write).
    inline const rPoint(dim)* getPoints() const { return pts; } //!< Accessor to the array of points.
    inline rPoint(dim) normal(int i) const { return normals[i]; } //!< Accessor to a normal (read).
    inline rPoint(dim)& normal(int i) { return normals[i]; } //!< Accessor to a normal (read/write).
    inline iPoint(dim) cell(int i) const { return cells[i]; } //!< Accessor to a cell (read).
    inline iPoint(dim)& cell(int i) { return cells[i]; } //!< Accessor to a cell (read/write).
    
    
#ifdef USE_VTK
    /*! \brief Clean the mesh.
        Remove points that are close (close means less than tolerance * maximum width of the mesh).\n
        Remove degenerate cells (= a cell made of one or two different points).\n
        Remove unused points.\n
        And recompute normals.
        \param tolerance Tolerance used to define what are two points close (default value: 1e-6).
        ONLY 3D !!!
    */
    void cleanMesh(float tolerance = 1e-6) {
       assert(dim==3);
       vtkPoints *points;
       vtkCellArray *polys;
       vtkPolyData *polyData;
       vtkCleanPolyData *cleaner;

       points=vtkPoints::New();
       points->SetNumberOfPoints(npts);
       for (int i=0; i<npts; i++)
           points->SetPoint(i,pts[i][0],pts[i][1],pts[i][2]);
       polys=vtkCellArray::New();
       polys->Allocate(ncls); // For speed
       for (int i=0; i<ncls; i++)
           polys->InsertNextCell(3,(int*)(cells+i));

       polyData=vtkPolyData::New();
       polyData->SetPoints(points);
       polyData->SetPolys(polys);
       polyData->Modified();

       cout<<"Number of points before cleaning mesh: "<<npts<<std::endl;;
       cout<<"Number of cells before cleaning mesh: "<<ncls<<std::endl;;

       cleaner=vtkCleanPolyData::New();
       cleaner->SetInput(polyData);
       cleaner->SetTolerance(tolerance);
       cleaner->ConvertLinesToPointsOn();
       cleaner->ConvertStripsToPolysOn();
       cleaner->CreateDefaultLocator();
       cleaner->Update();
    
       vtkPolyData* newpoly=cleaner->GetOutput();

       npts=newpoly->GetNumberOfPoints();
       ncls=newpoly->GetNumberOfCells();
       delete[] pts;
       delete[] cells;
       cells=new iPoint(3)[ncls];
    
       if (ncls)
       {
           int cpt=0;
    
           for (int i=0;i<ncls;i++)
           {
               if(newpoly->GetCell(i)->GetNumberOfPoints()==3)
               {
                   cells[cpt][0]=newpoly->GetCell(i)->GetPointId(0);
                   cells[cpt][1]=newpoly->GetCell(i)->GetPointId(1);
                   cells[cpt][2]=newpoly->GetCell(i)->GetPointId(2);
                   cpt++;
               }
           }
    
           if(cpt!=ncls)
           {
               iPoint(3) *newcells=new iPoint(3)[cpt];
               for(int i=0;i<cpt;i++) newcells[i]=cells[i];
               delete[] cells;
               cells=newcells;
               ncls=cpt;
           }
       }

       // Check if each point belongs to a cell
       int *sortTab=new int[npts]; for(int i=0;i<npts;i++) sortTab[i]=0;
       int *newIndex=new int[npts]; for(int i=0;i<npts;i++) newIndex[i]=0;
       for (int i=0;i<ncls;i++)
           for(int j=0;j<3;j++)
               sortTab[cells[i][j]]=1;

       int sum=-1;
       for(int i=0;i<npts;i++)
       {
           sum+=sortTab[i];
           newIndex[i]=sum;
       }
       
       // Update points
       npts=sum+1;
       pts=new rPoint(3)[npts];
       if (npts)
           for (int i=0;i<npts;i++) {
               if(sortTab[i]) {
                   pts[newIndex[i]][0]=newpoly->GetPoint(i)[0];
                   pts[newIndex[i]][1]=newpoly->GetPoint(i)[1];
                   pts[newIndex[i]][2]=newpoly->GetPoint(i)[2];
               }
           }
       
       // Update Cells
       for (int i=0;i<ncls;i++)
           for(int j=0;j<3;j++)
               cells[i][j]=newIndex[cells[i][j]];

       // Cleaning
       delete[] sortTab;
       delete[] newIndex;
       points->Delete();
       polys->Delete();
       polyData->Delete();
       cleaner->Delete();

       // Update normals
       if(withNormals) computeNormals();

       cout<<"Number of points after cleaning mesh: "<<npts<<std::endl;
       cout<<"Number of cells after cleaning mesh: "<<ncls<<std::endl;
       
    }
#endif

    /*! Recompute normals.
        \attention One normal for each point.
        ONLY 3D !!
    */
    void computeNormals() {
        assert(dim==3);
        if(withNormals) delete[] normals;
        normals=new rPoint(dim)[npts];
        withNormals=true;

        for (int i=0;i<npts;i++)
            normals[i]=rPoint(3)(0,0,0);

        for (int i=0;i<ncls;i++)
        {
            rPoint(3) n=(pts[cells[i][1]]-pts[cells[i][0]])^(pts[cells[i][2]]-pts[cells[i][0]]);
            n.normalize();
            for(int k=0;k<3;k++)
            {
                const rPoint(3) &A=pts[cells[i][k]];
                const rPoint(3) &B=pts[cells[i][(k+1)%3]];
                const rPoint(3) &C=pts[cells[i][(k+2)%3]];

                rPoint(3) BC=C-B;
                rPoint(3) BA=A-B;

                double angle=acos((BC*BA)/(BC.norme()*BA.norme()));
                normals[cells[i][(k+1)%3]]+=(angle*n);
            }
        }

        for (int i=0;i<npts;i++)
            normals[i].normalize();
    }

    /*! \brief Load a mesh.
        \attention Format recognized: VTK (from vtk), 3D (Nd), OFF, MESH (not implemented), TRI, MD (from robotvis ~xml, not implemented), DSGL (brainsuite).
        \param name filename.
        \param zoom zoom to apply to data (defaut value = 1).
        \param mft mesh file format (default value = NotSpecified means auto-recognize).
    */
    void load(const char* name,reel zoom=1,MeshFileType mft=NotSpecified) {
        char extension[128];
        getNameExtension(name,extension);
        if(!strcmp(extension,"vtk") || !strcmp(extension,"VTK") || mft==MFT_VTK) load_VTK(name,zoom);
        else if(!strcmp(extension,"3d") || !strcmp(extension,"3D")|| mft==MFT_3D) load_3D(name,zoom);
        else if(!strcmp(extension,"off") || !strcmp(extension,"OFF") || mft==MFT_OFF) load_OFF(name,zoom);
        else if(!strcmp(extension,"mesh") || !strcmp(extension,"MESH") || mft==MFT_MESH) load_MESH(name,zoom);
        else if(!strcmp(extension,"tri") || !strcmp(extension,"TRI") || mft==MFT_TRI) load_TRI(name,zoom);
        else if(!strcmp(extension,"md") || !strcmp(extension,"MD") || mft==MFT_MD) load_MD(name,zoom);
        else if(!strcmp(extension,"dsgl") || !strcmp(extension,"DSGL") || mft==MFT_DSGL) load_DSGL(name,zoom);
    }

    //! subfunction of load(const char*,reel,MeshFileType).
    void load_TRI(const char* name,reel zoom=1) {
        if (dim!=3)
        {
            std::cerr<<"    LsLib: mesh: TRI file format only in 3D"<<std::endl;
            return;
        }
        kill();
        count=new int(1);
        std::ifstream f(name);
        if(!f.is_open()) {std::cerr<<"Error opening file: "<<name<<std::endl; exit(1);}
        char ch; f>>ch;
        f>>npts;

        reel r;
        char myline[256];
        f.seekg( 0, std::ios_base::beg );
        f.getline(myline,256);
        f.getline(myline,256);
        int nread=0;
#if(USE_FLOATS)
        nread=sscanf(myline,"%f %f %f %f %f %f",&r,&r,&r,&r,&r,&r);
#else
        nread=sscanf(myline,"%lf %lf %lf %lf %lf %lf",&r,&r,&r,&r,&r,&r);
#endif
        withNormals= (nread==6);
        if(withNormals) normals = new rPoint(dim)[npts];
        f.seekg( 0, std::ios_base::beg );
        f.getline(myline,256);

        pts=new rPoint(dim)[npts];
        for (int i=0;i<npts;i++) {
            f>>pts[i];
            pts[i]=pts[i]*zoom;
            if(withNormals) f >> normals[i]; // Saute la normale (format Nicolas)
        }
        f>>ch;
        f>>ncls; f>>ncls; f>>ncls; // Ce numéro est répété tois fois.
        cells=new iPoint(dim)[ncls];
        for (int i=0;i<ncls;i++)
            f>>cells[i];
    }
    
    //! subfunction of load(const char*,reel,MeshFileType).
    void load_OFF(const char* name, reel zoom=1) {
        if (dim!=3) {
            std::cerr<<"    LsLib: mesh: OFF file format only in 3D"<<std::endl;
            return;
        }

        kill();
        count=new int(1);
        std::ifstream f(name);
        if(!f.is_open()) {std::cerr<<"Error opening file: "<<name<<std::endl; exit(1);}
        char poubelle[128]; int pb;
        f>>poubelle;        // put the "OFF" string
        f>>npts; f>>ncls; f>>pb;
        pts=new rPoint(dim)[npts];
        for (int i=0;i<npts;i++) {
            f>>pts[i];
            pts[i]=pts[i]*zoom;
        }
        cells=new iPoint(dim)[ncls];
        for (int i=0;i<ncls;i++) {
            f>>pb;        // put the "3"
            f>>cells[i];
        }
        withNormals=false;
    }

    /*! subfunction of load(const char*,reel,MeshFileType).
        \attention Only loads 3D meshes.
    */
    void load_VTK(const char* name, reel zoom=1) {
        if (dim!=3) {
            std::cerr << "    LsLib: mesh: VTK file format only in 3D for the moment ";
            std::cerr << "(others not implemented yet)" << std::endl;
            exit(1);
            return;
        }
#ifdef USE_VTK
        kill();
        count = new int(1);

        vtkPolyDataReader *reader = vtkPolyDataReader::New();
        reader->SetFileName(name);
        if (!reader->IsFilePolyData()) {
            std::cerr << "This is not a valid vtk poly data file" << std::endl;
            reader->Delete();
            exit(1);
        }
        reader->Update();
        vtkPolyData *vtkMesh = reader->GetOutput();
        reader->Delete();
        
        npts = vtkMesh->GetNumberOfPoints();

        pts=new rPoint(dim)[npts];
        
        for (int i=0;i<npts;i++)
        {
            pts[i][0]=vtkMesh->GetPoint(i)[0]*zoom;
            pts[i][1]=vtkMesh->GetPoint(i)[1]*zoom;
            pts[i][2]=vtkMesh->GetPoint(i)[2]*zoom;
        }
        ncls = vtkMesh->GetNumberOfCells();
        
        cells=new iPoint(dim)[ncls];
        vtkIdList *l;
        for (int i=0;i<ncls;i++) {
            if (vtkMesh->GetCellType(i) == VTK_TRIANGLE) {
                l = vtkMesh->GetCell(i)->GetPointIds();
                cells[i][0] = l->GetId(0);
                cells[i][1] = l->GetId(1);
                cells[i][2] = l->GetId(2);
            } else {
                std::cerr << "This is not a triangulation" << std::endl;
                exit(1);
            }
        }
        
        vtkMesh->Delete();
        l->Delete();

        withNormals=false; // TODO : read the normals if exists
#else
        std::cerr << "You have to compile OpenMEEG with VTK to read VTK files" << std::endl;
        exit(1);
#endif
    }

    /*! subfunction of load(const char*,reel,MeshFileType).
        \attention Loads Nd meshes (2D, 3D, 4D, etc.)
    */
    void load_3D (const char* name,reel zoom=1) {
        kill();
        count=new int(1);
        std::ifstream f(name);
        f>>npts;
        pts=new rPoint(dim)[npts];
        for (int i=0;i<npts;i++) {
            f>>pts[i];
            pts[i]=pts[i]*zoom;
            //rPoint(dim) foo; f >> foo; // Saute la normale (format Nicolas)
        }
        f>>ncls;
        cells=new iPoint(dim)[ncls];
        for (int i=0;i<ncls;i++)
            f>>cells[i];
        f >> withNormals;
        if (withNormals) {
            normals=new rPoint(dim)[npts];
            for (int i=0;i<npts;i++)
                f>>normals[i];
        }
        f.close();
    }

    /*! subfunction of load(const char*,reel,MeshFileType).
        \attention Not yet implemented.
    */
    void load_MESH(const char* name,reel zoom=1) {
        std::cerr<<"load_MESH not implemented yet (name: "<<name << " zoom:" << zoom <<")."<<std::endl;
        exit(1);
    }

    /*! subfunction of load(const char*,reel,MeshFileType).
        \attention Not yet implemented.
    */
    void load_MD(const char* name,reel zoom=1) {
        std::cerr<<"load_MD not implemented yet (name: "<<name << " zoom:" << zoom <<")."<<std::endl;
        exit(1);
    }

    //! subfunction of load(const char*,reel,MeshFileType).
    void load_DSGL(const char* name, reel zoom=1)
    {
        if (dim!=3)
        {
            std::cerr<<"    LsLib: mesh: DSGL file format only in 3D"<<std::endl;
            return;
        }
        if (zoom!=1) {
            std::cerr<<"load_DSGL not implemented yet with zoom != 1 (zoom = "<< zoom <<" )."<<std::endl;
        }
        kill();
        count=new int(1);

        FILE* fin=fopen(name,"rb");
        if(!fin) {
            std::cerr<<"Error opening file: "<<name<<std::endl;
            exit(1);
        }

        char *lgsd_string=(char*)"lgsd";
        char *dsgl_string=(char*)"dsgl";
        char chtab[5];
        fread(chtab,sizeof(char),4,fin);
        if(!strncmp(chtab,lgsd_string,4)) {std::cerr<<"dsgl big endian format not supported yet"<<std::endl; exit(1);}
        if(strncmp(chtab,dsgl_string,4)) {std::cerr<<"dsgl file format error"<<std::endl; exit(1);}
        int header_size=0;
        fread(&header_size,sizeof(int),1,fin);

        char version[8];
        if(header_size>44) {
            fread(version,sizeof(char),8,fin);
        }

        fread(&this->ncls,sizeof(int),1,fin);
        fread(&this->npts,sizeof(int),1,fin);

        float three_floats[3];
        float three_ints[3];
        fread(&three_ints,sizeof(int),1,fin);    //skip the nStripPoints
        fread(&three_floats,sizeof(float),3,fin);    //skip the resolutions
        if(header_size>44) {
            fread(&three_floats,sizeof(float),3,fin);    //skip the origin
        }
        fseek(fin,header_size,0);

        //read the faces
        this->cells=new iPoint(dim)[ncls];
        fread(cells,sizeof(int),3*ncls,fin);

        //read the vertices
        this->pts=new rPoint(dim)[npts];
        for(int i=0;i<npts;i++) {
            fread(&three_floats,sizeof(float),3,fin);
            this->pts[i]=rPoint(dim)(three_floats[0], three_floats[1], three_floats[2]);
        }
        fclose(fin);
        withNormals=false;  // no Normals in this file format
    }

    /*! \brief Save a mesh.
        \attention Format recognized: VTK (from vtk), 3D (Nd), OFF, MESH (not implemented), TRI, MD (from robotvis ~xml, not implemented).
        \param name filename.
        \param zoom zoom to apply to data (defaut value = 1).
        \param mft mesh file format (default value = NotSpecified means auto-recognize).
    */
    void save(const char* name,reel zoom=1,MeshFileType mft=NotSpecified) const {
        char extension[128];
        getNameExtension(name,extension);
        if(!strcmp(extension,"vtk") || !strcmp(extension,"VTK") || mft==MFT_VTK) save_VTK(name,zoom);
        else if(!strcmp(extension,"3d") || !strcmp(extension,"3D") || mft==MFT_3D) save_3D(name,zoom);
        else if(!strcmp(extension,"off") || !strcmp(extension,"OFF") || mft==MFT_OFF) save_OFF(name,zoom);
        else if(!strcmp(extension,"mesh") || !strcmp(extension,"MESH") || mft==MFT_MESH) save_MESH(name,zoom);
        else if(!strcmp(extension,"tri") || !strcmp(extension,"TRI") || mft==MFT_TRI) save_TRI(name,zoom);
        else if(!strcmp(extension,"md") || !strcmp(extension,"MD") || mft==MFT_MD) save_MD(name,zoom);
    }

    //! subfunction of save(const char*,reel,MeshFileType).
    void save_3D(const char* name,double zoom=1) const {
        std::ofstream f(name);
        f<<npts<<std::endl;
        for (int i=0;i<npts;i++)
            f<<pts[i]*zoom<<std::endl;
        f<<ncls<<std::endl;
        for (int i=0;i<ncls;i++)
            f<<cells[i]<<std::endl;
        f << withNormals << std::endl;
        if (withNormals)
            for (int i=0;i<npts;i++)
                f<<normals[i]<<std::endl;
        f.close();
    }

    /*! subfunction of save(const char*,reel,MeshFileType).
        \attention Save only 3D meshes.
    */
    void save_VTK ( const char* filename , reel zoom=1) const {
        if (dim!=3)
        {
            std::cerr << "    LsLib: mesh: VTK file format only in 3D for the moment ";
            std::cerr << "(others not implemented yet)" << std::endl;
            return;
        }
        std::ofstream os(filename);
        os<<"# vtk DataFile Version 2.0"<<std::endl;
        os<<"File "<<filename<<" generated by LsLib"<<std::endl;
        os<<"ASCII"<<std::endl;
        os<<"DATASET POLYDATA"<<std::endl;
        os<<"POINTS "<<npts<<" float"<<std::endl;
        for(int i=0;i<npts;i++) os<<pts[i]*zoom<<std::endl;
        os<<"POLYGONS "<<ncls<<" "<<ncls*4<<std::endl;
        for(int i=0;i<ncls;i++) os<<3<<" "<<cells[i]<<std::endl;
        os.close();
    }

    //! subfunction of save(const char*,reel,MeshFileType).
    void save_TRI ( const char* filename , reel zoom=1) const {
        if (dim!=3)
        {
            std::cerr<<"    LsLib: mesh: TRI file format only in 3D"<<std::endl;
            return;
        }
        bool save_withNormals=withNormals;
        if(!withNormals) ((mesh<3>*)this)->computeNormals();

        std::ofstream os(filename);
        os<<"- "<<this->npts<<std::endl;
        for(int i=0;i<npts;i++) {
            os<<this->pts[i]*zoom<<" ";
            os<<this->normals[i]<<std::endl;
        }
        os<<"- "<<ncls<<" "<<ncls<<" "<<ncls<<std::endl;
        for(int i=0;i<ncls;i++) os<<this->cells[i]<<std::endl;
        os.close();

        if(!save_withNormals) delete[] normals;
        *(bool*)(&withNormals)=save_withNormals;
    }

    /*! subfunction of save(const char*,reel,MeshFileType).
        \attention Save only 3D meshes.
    */
    void save_OFF ( const char* filename , reel zoom=1) const {
        if (dim!=3)
        {
            std::cerr<<"    LsLib: mesh: TRI file format only in 3D"<<std::endl;
            exit(1);
        }
        if (zoom!=1) {
            std::cerr<<"save_OFF not implemented yet with zoom!=1 (zoom = "<< zoom <<" )."<<std::endl;
            exit(1);
        }
        std::ofstream os(filename);
        os<<"OFF"<<std::endl;
        os<<npts<<" "<<ncls<<" "<<0<<std::endl;
        for(int i=0;i<npts;i++) os<<this->pts[i]<<std::endl;
        for(int i=0;i<ncls;i++) os<<(int)3<<" "<<this->cells[i]<<std::endl;
        os.close();
    }

    /*! subfunction of save(const char*,reel,MeshFileType).
        \attention Not implemented yet.
    */
    void save_MESH(const char* name,reel zoom=1) const {
        std::cerr<<"save_MESH not implemented yet (name: "<<name << " zoom:" << zoom <<")."<<std::endl;
        exit(1);
    }

    /*! subfunction of save(const char*,reel,MeshFileType).
        \attention Not implemented yet.
    */
    void save_MD(const char* name,reel zoom=1) const {
        std::cerr<<"save_MD not implemented yet (name: "<<name << " zoom:" << zoom <<")."<<std::endl;
        exit(1);
    }

    // Union
    /*! Create a mesh from two mesh.
        \attention No checking of common points and cells in both mesh.
    */
    mesh<dim> operator+(const mesh<dim>& M) const {
        if (!count)
            return M.clone();
        assert(withNormals==M.withNormals);
        mesh N(M.npts+npts,M.ncls+ncls,withNormals);
        for (int i=0;i<npts;i++)
            N.pts[i]=pts[i];
        for (int i=0;i<ncls;i++)
            N.cells[i]=cells[i];
        if (withNormals)
            for (int i=0;i<npts;i++)
                N.normals[i]=normals[i];
        for (int i=0;i<M.npts;i++)
            N.pts[npts+i]=M.pts[i];
        for (int i=0;i<M.ncls;i++)
            N.cells[ncls+i]=M.cells[i]+iPoint(dim)(npts);
        if (withNormals)
            for (int i=0;i<M.npts;i++)
                N.normals[npts+i]=M.normals[i];
        return N;
    }

    //! Remove normals.
    void killNormals() {
        if(withNormals)
        {
            withNormals=false;
            delete[] normals;
        }
    }

    /*! \param zoom Zoom factor.
        \return a clone of the meshes multiplied by the zoom factor.
    */
    mesh<dim> operator*(reel zoom) const {
        mesh<dim> N(clone());
        for (int i=0;i<N.npts;i++)
            N.pts[i]=N.pts[i]*zoom;
        return N;
    }

    /*! \param t translation factor
    */
    void translate(rPoint(dim) t) {
        for (int i=0;i<npts;i++)
            pts[i]+=t;
    }
    
    /*! \param s scale factor
    */
    void scale(reel s) {
        for (int i=0;i<npts;i++)
            pts[i]=pts[i]*s;
    }
    
    /*! \brief Return the list of cells linked to each point.
        \return an array of \c npts set of \c int .
    */
    intSet* pointCells() const {
        intSet* links=new intSet[npts];
        for (int i=0;i<ncls;i++)
            for (int d=0;d<dim;d++)
                links[cells[i][d]].insert(i);
        return links;
    }

    /*! \brief Return the list of cells linked to each cells (at least one common point).
        \return an array of \c ncls set of \c int .
    */
    intSet* cellCells() const {
        intSet* pc=pointCells();
        intSet* links=new intSet[ncls];
        for (int i=0;i<ncls;i++)
        {
            links[i]=pc[cells[i][0]];
            for (int d=1;d<dim;d++)
                links[i].insert( pc[cells[i][d]].begin() , pc[cells[i][d]].end() );
        }
        delete[] pc;
        return links;
    }

    /*! area of a cell
        \param i cell idx
    */
    reel size_cell(int i) const {
        std::cerr<<"size_cell not implemented yet for dim != 2 or 3"<<std::endl;
        return 0.;
    }

    /*! Sum of values on the Mesh
        \param values what to sum
    */
    template <typename T> T SumOnMesh(const T* values) const {
        T s(0.);
        for(int i=0;i<ncls;i++)
        {
            T t(0.);
            for(int k=0;k<dim;k++)
                t += values[cells[i][k]];

            s += t*(size_cell(i)/dim);
        }
        return(s);
    }

    /*! Average of values on the Mesh
        \param values what to compute average of
    */
    template <typename T> T AverageOnMesh(const T* values) const {
        T s(0.);
        reel length=0;
        for(int i=0;i<ncls;i++)
        {
            T t(0.);
            for(int k=0;k<dim;k++)
                t += values[cells[i][k]];
            reel delta=size_cell(i);
            length += delta;
            s += t*(delta/dim);
        }
        return(length?(s/length):T(0.));
    }

};

template<> reel mesh<2>::size_cell(int i) const
{
    return (point(cell(i)[0])-point(cell(i)[1])).norme();
}

template<> reel mesh<3>::size_cell(int i) const
{
    rPoint(3) P0=point(cell(i)[0]);
    rPoint(3) P1=point(cell(i)[1]);
    rPoint(3) P2=point(cell(i)[2]);
    return ((P0-P1)^(P2-P1)).norme()/2;
}

#endif
