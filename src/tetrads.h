/*
 * tetrads.h
 *
 *  Created on: Sep 12, 2019
 *      Author: bprather
 */

#ifndef TETRADS_H
#define TETRADS_H

#include "decs.h"

void coordinate_to_tetrad(REAL Ecov[NDIM][NDIM], REAL K[NDIM],
                      REAL K_tetrad[NDIM]);
void tetrad_to_coordinate(REAL Ecov[NDIM][NDIM], REAL K_tetrad[NDIM],
                      REAL K[NDIM]);
void set_Econ_from_trial(REAL Econ[4], int defdir, REAL trial[4]);
int check_handedness(REAL Econ[NDIM][NDIM], REAL Gcov[NDIM][NDIM], REAL *dot);
void project_out(REAL vcona[NDIM], REAL vconb[NDIM], REAL Gcov[4][4]);

#endif /* TETRADS_H */
