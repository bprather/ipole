/*
 * geometry.h
 *
 *  Created on: Sep 9, 2019
 *      Author: bprather
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "decs.h"

int gcon_func(REAL gcov[][NDIM], REAL gcon[][NDIM]);
REAL gdet_func(REAL gcov[][NDIM]);
void get_connection(REAL *X, REAL lconn[][NDIM][NDIM]);

void flip_index(REAL *ucon, REAL Gcov[NDIM][NDIM], REAL *ucov);
void flip_index_double(double *ucon, double Gcov[NDIM][NDIM], double *ucov);
// Old names aliased
inline void lower(REAL *ucon, REAL Gcov[NDIM][NDIM], REAL *ucov) {flip_index(ucon, Gcov, ucov);};
inline void raise(REAL *ucov, REAL Gcon[NDIM][NDIM], REAL *ucon) {flip_index(ucov, Gcon, ucon);};

void null_normalize(REAL Kcon[NDIM], REAL fnorm);
void normalize(REAL *vcon, REAL gcov[][NDIM]);


int invert_matrix(REAL Am[][NDIM], REAL Aminv[][NDIM]);
REAL theta_func(REAL X[NDIM]);
int levi_civita(int i, int j, int k, int l);

#endif /* GEOMETRY_H */
