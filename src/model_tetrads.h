/*
 * model_tetrads.h
 *
 *  Created on: Sep 12, 2019
 *      Author: bprather
 */

#ifndef MODEL_TETRADS_H
#define MODEL_TETRADS_H

#include "decs.h"

int make_camera_tetrad(REAL X[NDIM], REAL Econ[NDIM][NDIM],
                    REAL Ecov[NDIM][NDIM]);
int make_camera_tetrad_old(REAL X[NDIM], REAL Econ[NDIM][NDIM],
                        REAL Ecov[NDIM][NDIM]);
int make_plasma_tetrad(REAL Ucon[NDIM], REAL Kcon[NDIM], REAL Bcon[NDIM],
                    REAL Gcov[NDIM][NDIM], REAL Econ[NDIM][NDIM],
                    REAL Ecov[NDIM][NDIM]);

#endif /* MODEL_TETRADS_H */
