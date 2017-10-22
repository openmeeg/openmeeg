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

#include <assemble.h>
#include <danielsson.h>
#include <operators.h>
#include <sensors.h>
#include <om_utils.h>

namespace OpenMEEG {

    void assemble_ferguson(const Geometry& geo,Matrix& mat,const Matrix& pts);

    // EEG patches positions are reported line by line in the positions Matrix
    // mat is supposed to be filled with zeros
    // mat is the linear application which maps x (the unknown vector in symmetric system) -> v (potential at the electrodes)

    void
    assemble_Head2EEG(SparseMatrix& mat,const Geometry& geo,const Matrix& positions ) {

        mat = SparseMatrix(positions.nlin(),(geo.size()-geo.nb_current_barrier_triangles()));

        Vect3 current_position;
        Vect3 current_alphas;
        Triangle current_triangle;
        for ( unsigned i = 0; i < positions.nlin(); ++i) {
            for ( unsigned k = 0; k < 3; ++k) {
                current_position(k) = positions(i,k);
            }
            double dist;
            dist_point_geom(current_position,geo,current_alphas,current_triangle,dist);
            mat(i,current_triangle.s1().index()) = current_alphas(0);
            mat(i,current_triangle.s2().index()) = current_alphas(1);
            mat(i,current_triangle.s3().index()) = current_alphas(2);
        }
    }

    Head2EEGMat::Head2EEGMat(const Geometry& geo,const Sensors& electrodes) {
        assemble_Head2EEG(*this,geo,electrodes.getPositions());
    }

    // ECoG positions are reported line by line in the positions Matrix
    // mat is supposed to be filled with zeros
    // mat is the linear application which maps x (the unknown vector in symmetric system) -> v (potential at the ECoG electrodes)
    // difference with Head2EEG is that it interpolates the inner skull layer instead of the scalp layer. 

    void
    assemble_Head2ECoG(SparseMatrix& mat,const Geometry& geo,const Matrix& positions,const Interface& i) {
        mat = SparseMatrix(positions.nlin(),(geo.size()-geo.nb_current_barrier_triangles()));

        Vect3 current_position;
        Vect3 current_alphas;
        Triangle current_triangle;
        for (unsigned it = 0; it<positions.nlin(); ++it) {
            for (unsigned k = 0; k<3; ++k)
                current_position(k) = positions(it,k);
            dist_point_interface(current_position,i,current_alphas,current_triangle);
            mat(it,current_triangle.s1().index()) = current_alphas(0);
            mat(it,current_triangle.s2().index()) = current_alphas(1);
            mat(it,current_triangle.s3().index()) = current_alphas(2);
        }
    }

    Head2ECoGMat::Head2ECoGMat(const Geometry& geo,const Sensors& electrodes,const Interface& i) {
        assemble_Head2ECoG(*this,geo,electrodes.getPositions(),i);
    }

    Head2ECoGMat::Head2ECoGMat(const Geometry& geo,const Sensors& electrodes,const std::string& id) {
        assemble_Head2ECoG(*this,geo,electrodes.getPositions(),geo.interface(id));
    }

    // MEG patches positions are reported line by line in the positions Matrix (same for positions)
    // mat is supposed to be filled with zeros
    // mat is the linear application which maps x (the unknown vector in symmetric system) -> bFerguson (contrib to MEG response)

    void
    assemble_Head2MEG(Matrix& mat,const Geometry& geo,const Sensors& sensors) {

        Matrix positions = sensors.getPositions();
        Matrix orientations = sensors.getOrientations();
        const unsigned nbIntegrationPoints = sensors.getNumberOfPositions();
        unsigned p0_p1_size = (geo.size() - geo.nb_current_barrier_triangles());

        Matrix FergusonMat(3*nbIntegrationPoints,geo.nb_vertices());
        FergusonMat.set(0.0);

        assemble_ferguson(geo,FergusonMat,positions);

        mat = Matrix(nbIntegrationPoints,p0_p1_size);
        mat.set(0.0);

        for ( unsigned i = 0; i < nbIntegrationPoints; ++i) {
            PROGRESSBAR(i,nbIntegrationPoints);
            for ( Vertices::const_iterator vit = geo.vertex_begin(); vit != geo.vertex_end(); ++vit) {
                Vect3 fergusonField(FergusonMat(3*i,vit->index()),FergusonMat(3*i+1,vit->index()),FergusonMat(3*i+2,vit->index()));
                Vect3 normalizedDirection(orientations(i,0),orientations(i,1),orientations(i,2));
                normalizedDirection.normalize();
                mat(i,vit->index()) = fergusonField * normalizedDirection;
            }
        }
        mat = sensors.getWeightsMatrix() * mat; // Apply weights
    }

    Head2MEGMat::Head2MEGMat(const Geometry& geo,const Sensors& sensors) {
        assemble_Head2MEG(*this,geo,sensors);
    }

    // MEG patches positions are reported line by line in the positions Matrix (same for positions)
    // mat is supposed to be filled with zeros
    // mat is the linear application which maps x (the unknown vector in symmetric system) -> binf (contrib to MEG response)

    void
    assemble_SurfSource2MEG(Matrix& mat,const Mesh& sources_mesh,const Sensors& sensors) {

        Matrix positions = sensors.getPositions();
        Matrix orientations = sensors.getOrientations();
        const unsigned nsquids = positions.nlin();

        mat = Matrix(nsquids,sources_mesh.nb_vertices());
        mat.set(0.0);

        for ( unsigned i = 0; i < nsquids; ++i) {
            PROGRESSBAR(i,nsquids);
            Vect3 p(positions(i,0),positions(i,1),positions(i,2));
            Matrix FergusonMat(3,mat.ncol());
            FergusonMat.set(0.0);
            operatorFerguson(p,sources_mesh,FergusonMat,0,1.);
            for ( unsigned j = 0; j < mat.ncol(); ++j) {
                Vect3 fergusonField(FergusonMat(0,j),FergusonMat(1,j),FergusonMat(2,j));
                Vect3 normalizedDirection(orientations(i,0),orientations(i,1),orientations(i,2));
                normalizedDirection.normalize();
                mat(i,j) = fergusonField * normalizedDirection;
            }
        }

        mat = sensors.getWeightsMatrix() * mat; // Apply weights
    }

    SurfSource2MEGMat::SurfSource2MEGMat(const Mesh& sources_mesh,const Sensors& sensors) {
        assemble_SurfSource2MEG(*this,sources_mesh,sensors);
    }

    // Creates the DipSource2MEG Matrix with unconstrained orientations for the sources.
    // MEG patches positions are reported line by line in the positions Matrix (same for positions)
    // mat is supposed to be filled with zeros
    // sources is the name of a file containing the description of the sources - one dipole per line: x1 x2 x3 n1 n2 n3,x being the position and n the orientation.

    void
    assemble_DipSource2MEG(Matrix& mat,const Matrix& dipoles,const Sensors& sensors) {

        Matrix positions = sensors.getPositions();
        Matrix orientations = sensors.getOrientations();

        if ( dipoles.ncol() != 6) {
            std::cerr << "Dipoles File Format Error" << std::endl;
            exit(1);
        }

        // this Matrix will contain the field generated at the location of the i-th squid by the j-th source
        mat = Matrix(positions.nlin(),dipoles.nlin());

        // the following routine is the equivalent of operatorFerguson for pointlike dipoles.
        for ( unsigned i = 0; i < mat.nlin(); ++i) {
            for ( unsigned j = 0; j < mat.ncol(); ++j) {
                Vect3 r(dipoles(j,0),dipoles(j,1),dipoles(j,2));
                Vect3 q(dipoles(j,3),dipoles(j,4),dipoles(j,5));
                Vect3 diff(positions(i,0),positions(i,1),positions(i,2));
                diff -= r;
                double norm_diff = diff.norm();
                Vect3 fergusonField = q ^ diff / (norm_diff * norm_diff * norm_diff);
                Vect3 normalizedDirection(orientations(i,0),orientations(i,1),orientations(i,2));
                normalizedDirection.normalize();
                mat(i,j) = fergusonField * normalizedDirection * MU0 / (4.0 * M_PI);
            }
        }

        mat = sensors.getWeightsMatrix() * mat; // Apply weights
    }

    DipSource2MEGMat::DipSource2MEGMat(const Matrix& dipoles,const Sensors& sensors) {
        assemble_DipSource2MEG(*this,dipoles,sensors);
    }

}
