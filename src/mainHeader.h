#ifndef H_mainHeader
#define H_mainHeader


#include <vector>
#include <vecteur.h>
#include <matrice.h>
#include <symmatrice.h>

#include "geometry.h"

void assemble_matrice(const Geometry &geo,symmatrice &mat,const int);
void assemble_RHSmatrix(const Geometry &geo,const Mesh &sources,matrice &mat,const int);
void assemble_RHS2matrix(const Geometry &geo,const Mesh &sources, matrice &mat,const int);
void assemble_RHSvector(const Geometry &geo, std::vector<Vect3> Rs, std::vector<Vect3> Qs, vecteur &rhs,const int);
void assemble_RHS_dipoles_matrice(const Geometry &geo, std::vector<Vect3> Rs, std::vector<Vect3> Qs, matrice &rhs,const int);
void assemble_RHS_dipoles_matrice_grad(const Geometry &geo, std::vector<Vect3> Rs, std::vector<Vect3> Qs, matrice &rhs,const int);
void assemble_ferguson(const Geometry &geo, matrice &mat, const Vect3* pts,const int n);
void assemble_xToEEGresponse(const Geometry &geo, matrice &mat, const matrice &positions );
void assemble_xToMEGresponseContrib(const Geometry &geo, matrice &mat, const matrice &positions, const matrice &orientations );
void assemble_sToMEGresponseContrib(const Mesh &sources_mesh, matrice &mat, const matrice &positions, const matrice &orientations );
void deflat(genericMatrix &M, int start, int end, double coef);

void assemble_sToMEGresponseContrib_point(matrice&dipoles, matrice &mat, const matrice &positions, const matrice &orientations );

#endif
