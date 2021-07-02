#ifndef SYMPHONY_POWER_LAW_H_
#define SYMPHONY_POWER_LAW_H_
#include "params.h"
//#include "distribution_function_common_routines.h"

REAL power_law_to_be_normalized(REAL gamma, void * paramsInput);
REAL power_law_f(REAL gamma, struct parameters * params); 
REAL differential_of_power_law(REAL gamma, struct parameters * params);

REAL power_law_I(struct parameters * params);
REAL power_law_Q(struct parameters * params);
REAL power_law_V(struct parameters * params);

REAL power_law_I_abs(struct parameters * params);
REAL power_law_Q_abs(struct parameters * params);
REAL power_law_V_abs(struct parameters * params);

#endif /* SYMPHONY_POWER_LAW_H_ */
