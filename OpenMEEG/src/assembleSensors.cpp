// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <assemble.h>
#include <danielsson.h>
#include <operators.h>
#include <sensors.h>

#include <constants.h>
#include <sparse_matrix.h>
#include <OMExceptions.H>

namespace OpenMEEG {

    void assemble_ferguson(const Geometry& geo,Matrix& mat,const Matrix& pts);

    // EEG patches positions are reported line by line in the positions Matrix
    // mat is supposed to be filled with zeros
    // mat is the linear application which maps x (the unknown vector in symmetric system) -> v (potential at the electrodes)

    SparseMatrix Head2EEGMat(const Geometry& geo,const Sensors& electrodes) {
        const Matrix& positions = electrodes.getPositions();
        SparseMatrix mat(positions.nlin(),(geo.nb_parameters()-geo.nb_current_barrier_triangles()));

        for (unsigned i=0;i<positions.nlin();++i) {
            const Vect3 current_position(positions(i,0),positions(i,1),positions(i,2));
            Vect3 current_alphas;
            const auto& res = dist_point_geom(current_position,geo,current_alphas);
            const Triangle& current_triangle = std::get<1>(res);
            for (unsigned j=0;j<3;++j)
                mat(i,current_triangle.vertex(j).index()) = current_alphas(j);
        }

        return mat;
    }

    // ECoG positions are reported line by line in the positions Matrix
    // mat is supposed to be filled with zeros
    // mat is the linear application which maps x (the unknown vector in symmetric system) -> v (potential at the ECoG electrodes)
    // difference with Head2EEG is that it interpolates the inner skull layer instead of the scalp layer.

    SparseMatrix Head2ECoGMat(const Geometry& geo,const Sensors& electrodes,const Interface& i) {

        const Matrix& positions = electrodes.getPositions();
        SparseMatrix mat(positions.nlin(),(geo.nb_parameters()-geo.nb_current_barrier_triangles()));

        for (unsigned it=0;it<positions.nlin();++it) {
            Vect3 current_position;
            for (unsigned k=0;k<3;++k)
                current_position(k) = positions(it,k);
            Vect3 current_alphas;
            const auto& res = dist_point_interface(current_position,i,current_alphas);
            const Triangle& current_triangle = std::get<1>(res);
            for (unsigned j=0;j<3;++j)
                mat(it,current_triangle.vertex(j).index()) = current_alphas(j);
        }

        return mat;
    }

    // MEG patches positions are reported line by line in the positions Matrix (same for positions)
    // mat is supposed to be filled with zeros
    // mat is the linear application which maps x (the unknown vector in symmetric system) -> bFerguson (contrib to MEG response)

    Matrix Head2MEGMat(const Geometry& geo,const Sensors& sensors) {

        const Matrix& positions    = sensors.getPositions();
        const Matrix& orientations = sensors.getOrientations();
        const unsigned nbIntegrationPoints = sensors.getNumberOfPositions();
        unsigned p0_p1_size = geo.nb_parameters()-geo.nb_current_barrier_triangles();

        Matrix FergusonMat(3*nbIntegrationPoints,geo.vertices().size());
        FergusonMat.set(0.0);

        assemble_ferguson(geo,FergusonMat,positions);

        Matrix mat(nbIntegrationPoints,p0_p1_size);
        mat.set(0.0);

        ProgressBar pb(nbIntegrationPoints);
        for (unsigned i=0; i<nbIntegrationPoints; ++i,++pb) {
            for (const auto& vertex : geo.vertices()) {
                const Vect3 fergusonField(FergusonMat(3*i,vertex.index()),
                                          FergusonMat(3*i+1,vertex.index()),
                                          FergusonMat(3*i+2,vertex.index()));
                const Vect3 direction(orientations(i,0),orientations(i,1),orientations(i,2));
                mat(i,vertex.index()) = dotprod(fergusonField,direction)/direction.norm();
            }
        }

        return sensors.getWeightsMatrix()*mat; // Apply weights
    }

    // MEG patches positions are reported line by line in the positions Matrix (same for positions)
    // mat is supposed to be filled with zeros
    // mat is the linear application which maps x (the unknown vector in symmetric system) -> binf (contrib to MEG response)

    Matrix SurfSource2MEGMat(const Mesh& sources_mesh,const Sensors& sensors) {

        const Matrix& positions    = sensors.getPositions();
        const Matrix& orientations = sensors.getOrientations();
        const unsigned nsquids = positions.nlin();

        Matrix mat(nsquids,sources_mesh.vertices().size());
        mat.set(0.0);

        ProgressBar pb(nsquids);
        for (unsigned i=0; i<nsquids; ++i,++pb) {
            Vect3 p(positions(i,0),positions(i,1),positions(i,2));
            Matrix FergusonMat(3,mat.ncol());
            FergusonMat.set(0.0);
            operatorFerguson(p,sources_mesh,FergusonMat,0,1.);
            for (unsigned j=0;j<mat.ncol();++j) {
                const Vect3 fergusonField(FergusonMat(0,j),FergusonMat(1,j),FergusonMat(2,j));
                const Vect3 direction(orientations(i,0),orientations(i,1),orientations(i,2));
                mat(i,j) = dotprod(fergusonField,direction)/direction.norm();
            }
        }

        return sensors.getWeightsMatrix()*mat; // Apply weights
    }

    // Creates the DipSource2MEG Matrix with unconstrained orientations for the sources.
    // MEG patches positions are reported line by line in the positions Matrix (same for positions)
    // mat is supposed to be filled with zeros
    // sources is the name of a file containing the description of the sources - one dipole
    // per line: x1 x2 x3 n1 n2 n3,x being the position and n the orientation.

    Matrix DipSource2MEGMat(const Matrix& dipoles,const Sensors& sensors) {

        const Matrix& positions    = sensors.getPositions();
        const Matrix& orientations = sensors.getOrientations();

        if (dipoles.ncol()!=6)
            throw OpenMEEG::DipoleError("expecting 6 columns.");

        // This Matrix will contain the field generated at the location of the i-th squid by the j-th source

        Matrix mat(positions.nlin(),dipoles.nlin());

        // The following routine is the equivalent of operatorFerguson for point-like dipoles.

        for (unsigned i=0;i<mat.nlin();++i)
            for (unsigned j=0;j<mat.ncol();++j) {
                const Vect3 r(dipoles(j,0),dipoles(j,1),dipoles(j,2));
                const Vect3 q(dipoles(j,3),dipoles(j,4),dipoles(j,5));
                const Vect3& diff = Vect3(positions(i,0),positions(i,1),positions(i,2))-r;
                const double norm_diff = diff.norm();
                const Vect3 fergusonField = q ^ diff / (norm_diff * norm_diff * norm_diff);
                const Vect3 direction(orientations(i,0),orientations(i,1),orientations(i,2));
                mat(i,j) = dotprod(fergusonField,direction)*MagFactor/direction.norm();
            }

        return sensors.getWeightsMatrix()*mat; // Apply weights
    }
}
