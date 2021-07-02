/*
 * radiation.h
 *
 *  Created on: Sep 12, 2019
 *      Author: bprather
 */

#ifndef RADIATION_H
#define RADIATION_H

#include "decs.h"

/* radiation */
REAL Bnu_inv(REAL nu, REAL Thetae);
REAL jnu_inv(REAL nu, REAL Thetae, REAL Ne, REAL B, REAL theta);
REAL get_fluid_nu(REAL Kcon[NDIM], REAL Ucov[NDIM]);
REAL get_bk_angle(REAL X[NDIM], REAL Kcon[NDIM], REAL Ucov[NDIM],
              REAL Bcon[NDIM], REAL Bcov[NDIM]);

#endif /* RADIATION_H */
