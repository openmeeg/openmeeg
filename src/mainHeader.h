#ifndef H_mainHeader
#define H_mainHeader


#include <vector>
#include <vecteur.h>
#include <matrice.h>
#include <symmatrice.h>

#include "geometry.h"

void assemble_matrice( Geometry &geo, symmatrice &mat);
void assemble_RHSmatrix(Geometry &geo, Mesh &sources, matrice &mat);
void assemble_RHS2matrix(Geometry &geo, Mesh &sources, matrice &mat);
void assemble_RHSvector( Geometry &geo, std::vector<Vect3> Rs, std::vector<Vect3> Qs, vecteur &rhs);
void assemble_RHS_dipoles_matrice( Geometry &geo, std::vector<Vect3> Rs, std::vector<Vect3> Qs, matrice &rhs);
void assemble_RHS_dipoles_matrice_grad( Geometry &geo, std::vector<Vect3> Rs, std::vector<Vect3> Qs, matrice &rhs);
void assemble_ferguson( Geometry &geo, matrice &mat, const Vect3* pts, int n);
void assemble_xToEEGresponse( Geometry &geo, matrice &mat, const matrice &positions );
void assemble_xToMEGresponseContrib( Geometry &geo, matrice &mat, const matrice &positions, const matrice &orientations );
void assemble_sToMEGresponseContrib( Mesh &sources_mesh, matrice &mat, const matrice &positions, const matrice &orientations );
void deflat(genericMatrix &M, int start, int end, double coef);

void assemble_sToMEGresponseContrib_point(matrice&dipoles, matrice &mat, const matrice &positions, const matrice &orientations );

#endif
