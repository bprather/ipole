#ifndef GRID_H
#define GRID_H

#include "decs.h"

extern int N1, N2, N3;

REAL gdet_zone(int i, int j, int k);

void ijktoX(int i, int j, int k, REAL X[NDIM]);
void Xtoijk(REAL X[NDIM], int *i, int *j, int *k, REAL del[NDIM]);

int X_in_domain(REAL X[NDIM]);

REAL interp_scalar(REAL X[NDIM], double ***var);

#endif // GRID_H
