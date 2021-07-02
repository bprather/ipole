
#ifndef MODEL_PARAMS_H
#define MODEL_PARAMS_H

#include "decs.h"

// There is really only one problem where we should record Stokes parameters per step
// it makes no sense in a relativistic context
#define INTEGRATOR_TEST (1)

// Model parameters (TODO eliminate if we can these are unused)
extern REAL rmax_geo;
extern REAL model_dl;

void record_stokes_parameters(REAL SI, REAL SQ, REAL SU, REAL SV, REAL lam);

#endif /* MODEL_PARAMS_H */
