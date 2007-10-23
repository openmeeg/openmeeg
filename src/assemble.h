#ifndef _ASSEMBLE_H_
#define _ASSEMBLE_H_

#include <vector>

#include "vecteur.h"
#include "matrice.h"
#include "symmatrice.h"
#include "geometry.h"

void assemble_LHS(const Geometry &geo,symmatrice &mat,const int);
void assemble_RHS(const Geometry &geo,const Mesh &sources,matrice &mat,const int);
void assemble_RHS2(const Geometry &geo,const Mesh &sources, matrice &mat,const int);
void assemble_RHS_dipoles(const Geometry &geo, std::vector<Vect3> Rs, std::vector<Vect3> Qs, matrice &rhs,const int);
void assemble_RHS_dipoles_grad(const Geometry &geo, std::vector<Vect3> Rs, std::vector<Vect3> Qs, matrice &rhs,const int);
void assemble_ferguson(const Geometry &geo, matrice &mat, const Vect3* pts,const int n);
void assemble_vToEEG(const Geometry &geo, matrice &mat, const matrice &positions );
void assemble_vToMEG(const Geometry &geo, matrice &mat, const matrice &positions, const matrice &orientations );
void assemble_sToMEG(const Mesh &sources_mesh, matrice &mat, const matrice &positions, const matrice &orientations );
void deflat(genericMatrix &M, int start, int end, double coef);

void assemble_sToMEG_point(matrice&dipoles, matrice &mat, const matrice &positions, const matrice &orientations );

#endif /* _ASSEMBLE_H_ */
