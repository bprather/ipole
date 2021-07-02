#ifndef COORDINATES_H
#define COORDINATES_H

#include "decs.h"

// Poor man's enum of coordinate systems
// Modified Kerr-Schild coordinates.  Gammie '03.
#define METRIC_MKS 0
// New-style BHAC MKS.  Just like MKS but with X2 in [0,pi] and hslope->(1-hslope)
#define METRIC_BHACMKS 1
// "Funky" MKS coordinates as defined on IL wiki, see
// https://github.com/AFD-Illinois/docs/wiki/Coordinates
#define METRIC_FMKS 2
// MKS3 coordinates from koral-light
#define METRIC_MKS3 3
// Spherical coordinates in Minkowski space
#define METRIC_MINKOWSKI 4

// Coordinate parameters.  See 
extern int use_eKS_internal;
extern int metric;
extern REAL a, hslope; // mks
extern REAL poly_norm, poly_xt, poly_alpha, mks_smooth; // fmks
extern REAL mks3R0, mks3H0, mks3MY1, mks3MY2, mks3MP0; // mks3
extern REAL startx[NDIM], stopx[NDIM], dx[NDIM];
extern REAL R0, Rin, Rout, Rh;

void bl_coord(REAL *X, REAL *r, REAL *th);
void bl_to_ks(REAL X[NDIM], REAL ucon_bl[NDIM], REAL ucon_ks[NDIM]);
void ks_to_bl(REAL X[NDIM], REAL ucon_ks[NDIM], REAL ucon_bl[NDIM]);
void gcov_func(REAL *X, REAL gcov[][NDIM]);
// TODO privatize these, why are they needed in models?
void gcov_ks(REAL r, REAL th, REAL gcov[NDIM][NDIM]);
void gcov_bl(REAL r, REAL th, REAL gcov[NDIM][NDIM]);

// Internal
void set_dxdX(REAL X[NDIM], REAL dxdX[NDIM][NDIM]);
void set_dXdx(REAL X[NDIM], REAL dXdx[NDIM][NDIM]);
void vec_to_ks(REAL X[NDIM], REAL v_nat[NDIM], REAL v_ks[NDIM]);
void vec_from_ks(REAL X[NDIM], REAL v_ks[NDIM], REAL v_nat[NDIM]);

// Translation to native coords
void native_coord(REAL r, REAL th, REAL phi, REAL X[NDIM]);
REAL root_find(REAL X[NDIM]);

#endif // COORDINATES_H
