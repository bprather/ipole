
#ifndef MODEL_PARAMS_H
#define MODEL_PARAMS_H

#include "decs.h"

#define THIN_DISK (1)

// Model parameters
extern REAL rmax_geo;

// New definitions needed for problem defined with boundary condition
// Uses are #ifdef'd off on the THIN_DISK flag above
int thindisk_region(REAL Xi[NDIM], REAL Xf[NDIM]);
void get_model_i(REAL X[NDIM], REAL K[NDIM], REAL *SI);
void get_model_stokes(REAL X[NDIM], REAL K[NDIM], REAL *SI, REAL *SQ, REAL *SU, REAL *SV);


#endif /* MODEL_PARAMS_H */
