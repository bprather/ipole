#ifndef MODEL_H
#define MODEL_H

#include "model_params.h"
#include "decs.h"
#include "par.h"
#include "hdf5_utils.h"

extern REAL M_unit;
extern REAL L_unit;
extern REAL T_unit;
extern REAL RHO_unit;
extern REAL U_unit;
extern REAL B_unit;
extern REAL Te_unit;

void try_set_model_parameter(const char *word, const char *value);
void init_model(REAL *tA, REAL *tB);

int radiating_region(REAL X[NDIM]);

double get_model_thetae(REAL X[NDIM]);
double get_model_b(REAL X[NDIM]);
double get_model_ne(REAL X[NDIM]);

// For exotic or custom distributions
void get_model_powerlaw_vals(REAL X[NDIM], REAL *p, REAL *n,
                             REAL *gamma_min, REAL *gamma_max, REAL *gamma_cut);
void get_model_jar(REAL X[NDIM], REAL Kcon[NDIM],
    REAL *jI, REAL *jQ, REAL *jU, REAL *jV,
    REAL *aI, REAL *aQ, REAL *aU, REAL *aV,
    REAL *rQ, REAL *rU, REAL *rV);
void get_model_jk(REAL X[NDIM], REAL Kcon[NDIM], REAL *jnuinv, REAL *knuinv);

void get_model_fourv(REAL X[NDIM], REAL Kcon[NDIM],
                     REAL Ucon[NDIM], REAL Ucov[NDIM],
                     REAL Bcon[NDIM], REAL Bcov[NDIM]);

void output_hdf5();

// Optional function to be able to trace primitives along a geodesic
void get_model_primitives(REAL X[NDIM], REAL *p);

#endif // MODEL_H
