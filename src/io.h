/*
 * io.h
 *
 *  Created on: Sep 10, 2019
 *      Author: bprather
 */

#ifndef IO_H
#define IO_H

#include "decs.h"
#include "par.h"

void dump(REAL image[], REAL imageS[], REAL taus[],
          const char *fname, REAL scale, REAL cam[NDIM],
          REAL fovx, REAL fovy, size_t nx, size_t ny, Params *params, int nopol);
void dump_var_along(int i, int j, int nstep, struct of_traj *traj, int nx, int ny,
                    REAL scale, REAL cam[NDIM], REAL fovx, REAL fovy, Params *params);

#endif // IO_H
