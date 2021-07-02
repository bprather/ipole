/*
 * model_geodesics.h
 *
 *  Created on: Sep 11, 2019
 *      Author: bprather
 */

#ifndef SRC_MODEL_GEODESICS_H_
#define SRC_MODEL_GEODESICS_H_

#include "decs.h"
#include "par.h"

int trace_geodesic(REAL Xi[NDIM], REAL Kconi[NDIM], struct of_traj *traj, REAL eps, int step_max);
void init_XK(long int i, long int j, int nx, int ny, REAL Xcam[NDIM],
             Params params, REAL fovx, REAL fovy,
             REAL X[NDIM], REAL Kcon[NDIM]);

// Internal utilities still used for slow light
int stop_backward_integration(REAL X[NDIM], REAL Xhalf[NDIM], REAL Kcon[NDIM]);
REAL stepsize(REAL X[NDIM], REAL K[NDIM], REAL eps);

#endif /* SRC_MODEL_GEODESICS_H_ */
