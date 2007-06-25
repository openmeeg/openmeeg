#ifndef H_mainHeader
#define H_mainHeader


#include <vector>
#include <vecteur.h>
#include <matrice.h>
#include <symmatrice.h>

#include "geometry.h"

void assemble_matrice( geometry &geo, symmatrice &mat);
void assemble_RHSmatrix(geometry &geo, mesh &sources, matrice &mat);
void assemble_RHS2matrix(geometry &geo, mesh &sources, matrice &mat);
void assemble_RHSvector( geometry &geo, std::vector<vect3> Rs, std::vector<vect3> Qs, vecteur &rhs);
void assemble_RHS_dipoles_matrice( geometry &geo, std::vector<vect3> Rs, std::vector<vect3> Qs, matrice &rhs);
void assemble_RHS_dipoles_matrice_grad( geometry &geo, std::vector<vect3> Rs, std::vector<vect3> Qs, matrice &rhs);
void assemble_ferguson( geometry &geo, matrice &mat, const vect3* pts, int n);
void assemble_xToEEGresponse( geometry &geo, matrice &mat, const matrice &positions );
void assemble_xToMEGresponseContrib( geometry &geo, matrice &mat, const matrice &positions, const matrice &orientations );
void assemble_sToMEGresponseContrib( mesh &sources_mesh, matrice &mat, const matrice &positions, const matrice &orientations );
void deflat(genericMatrix &M, int start, int end, double coef);

void assemble_sToMEGresponseContrib_point(matrice&dipoles, matrice &mat, const matrice &positions, const matrice &orientations );

#endif
