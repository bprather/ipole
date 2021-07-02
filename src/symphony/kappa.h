#ifndef SYMPHONY_KAPPA_H_
#define SYMPHONY_KAPPA_H_
#include "params.h"
//#include "distribution_function_common_routines.h"
#include "gsl/gsl_sf_hyperg.h"

REAL kappa_to_be_normalized(REAL gamma, void * paramsInput);
REAL kappa_f(REAL gamma, struct parameters * params);
REAL differential_of_kappa(REAL gamma, struct parameters * params);

REAL kappa_I(struct parameters * params);
REAL kappa_Q(struct parameters * params);
REAL kappa_V(struct parameters * params);

REAL kappa_I_abs(struct parameters * params);
REAL kappa_Q_abs(struct parameters * params);
REAL kappa_V_abs(struct parameters * params);

#endif /* SYMPHONY_KAPPA_H_ */
