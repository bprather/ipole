/*
 * ipolarray.h
 *
 *  Created on: Sep 12, 2019
 *      Author: bprather
 */

#ifndef IPOLARRAY_H
#define IPOLARRAY_H

#include "decs.h"
#include "par.h"

#include <complex.h>

// Top-level functions for solving emission
int integrate_emission(struct of_traj *traj, int nstep,
                    REAL *Intensity, REAL *Tau, REAL *tauF,
                    _Complex REAL N_coord[NDIM][NDIM], Params *params);

// Needed for slow light.  TODO extend above to use instead
int evolve_N(REAL Xi[NDIM],REAL Kconi[NDIM],
    REAL Xhalf[NDIM], REAL Kconhalf[NDIM],
    REAL Xf[NDIM],REAL Kconf[NDIM],
    REAL dlam,
    _Complex REAL N_coord[NDIM][NDIM],
    REAL *tauF, Params *params);
REAL approximate_solve (REAL Ii, REAL ji, REAL ki, REAL jf, REAL kf,
                   REAL dl, REAL *tau);

void project_N(REAL X[NDIM],REAL Kcon[NDIM],
    _Complex REAL Ncon[NDIM][NDIM],
    REAL *Stokes_I, REAL *Stokes_Q,REAL *Stokes_U,REAL *Stokes_V, REAL rotcam);

#endif /* IPOLARRAY_H */
