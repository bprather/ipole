#ifndef SYMPHONY_MAXWELL_JUETTNER_H_
#define SYMPHONY_MAXWELL_JUETTNER_H_
#include "params.h"
//#include "distribution_function_common_routines.h"
#include <gsl/gsl_sf_bessel.h>

REAL maxwell_juettner_f(REAL gamma, struct parameters * params);
REAL differential_of_maxwell_juettner(REAL gamma, struct parameters * params);
REAL maxwell_juettner_n_peak(struct parameters * params);

/* Fits */

/* Emissivities */
REAL maxwell_juettner_I(struct parameters * params);
REAL maxwell_juettner_Q(struct parameters * params);
REAL maxwell_juettner_V(struct parameters * params);

REAL planck_func(struct parameters * params);

/* Absorptivities */
REAL maxwell_juettner_I_abs(struct parameters * params);
REAL maxwell_juettner_Q_abs(struct parameters * params);
REAL maxwell_juettner_V_abs(struct parameters * params);

/* Faraday rotation/conversion coefficients */
REAL maxwell_juettner_rho_Q(struct parameters * params);
REAL maxwell_juettner_rho_V(struct parameters * params);

#endif /* SYMPHONY_MAXWELL_JUETTNER_H_ */
