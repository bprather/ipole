#ifndef DECS_H
#define DECS_H

#include "constants.h"
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>

// Strings and string tools
#define VERSION_STRING "ipole-beta-1.2"
#define xstr(s) str(s)
#define str(s) #s
#define STRLEN (2048)

// Dimensions
#define NDIM	4
#define NIMG (5) // Stokes parameters plus faraday depth

// "functions"
#define sign(x) (((x) < 0) ? -1 : ((x) > 0))
#define delta(x, y) ((x) == (y))

// MAX and MIN were here but are both unsafe & slow. Use fmax/fmin

#define MULOOP for(int mu=0;mu<NDIM;mu++)
#define MUNULOOP for(int mu=0;mu<NDIM;mu++) for(int nu=0;nu<NDIM;nu++)

#define H5T_DOUBLE H5T_IEEE_F64LE
// accurate...
//#define REAL _Float128
#define REAL long double
// ...or fast
//#define REAL double

struct of_traj {
  REAL dl;
  REAL X[NDIM];
  REAL Kcon[NDIM];
  REAL Xhalf[NDIM];
  REAL Kconhalf[NDIM];
};

// used for slow light
struct of_image {
  int nstep;
  REAL intensity;
  REAL tau;
  REAL tauF;
  _Complex REAL N_coord[NDIM][NDIM];
};

#endif // DECS_H
