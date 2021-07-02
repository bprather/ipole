/*
 * model_radiation.h
 *
 *  Created on: Sep 12, 2019
 *      Author: bprather
 */

#ifndef MODEL_RADIATION_H
#define MODEL_RADIATION_H

#include "decs.h"
#include "par.h"

/* transfer coefficients in tetrad frame */
void jar_calc(REAL X[NDIM], REAL Kcon[NDIM], REAL *jI, REAL *jQ,
              REAL *jU, REAL *jV, REAL *aI, REAL *aQ, REAL *aU,
              REAL *aV, REAL *rQ, REAL *rU, REAL *rV, Params *params);

REAL jnu_synch(REAL nu, REAL Ne, REAL Thetae, REAL B, REAL theta);
void get_jkinv(REAL X[NDIM], REAL Kcon[NDIM], REAL *jnuinv, REAL *knuinv, Params *params);

#endif /* MODEL_RADIATION_H */
