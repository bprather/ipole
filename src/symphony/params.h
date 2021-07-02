#ifndef SYMPHONY_PARAMS_H_
#define SYMPHONY_PARAMS_H_

#include "decs.h"

#include <math.h>

struct parameters
{
  /*parameters of calculation*/
  /*we use Gaussian CGS units*/
  REAL pi;        
  REAL mass_electron;
  REAL plancks_constant;
  REAL speed_light;
  REAL electron_charge;
  REAL n_max;
  int    C;       
  /*Keys for the distributions*/
  int    MAXWELL_JUETTNER;
  int    POWER_LAW;
  int    KAPPA_DIST;
  /*Keys for the polarization modes*/
  int    STOKES_I;
  int    STOKES_Q;
  int    STOKES_U;
  int    STOKES_V;
  /*Keys for the mode: absorptivity or emissivity*/
  int    ABSORPTIVITY;
  int    EMISSIVITY;
  /*Keys for alpha_nu computation method*/
  int SYMPHONY_METHOD;
  int SUSCEPT_METHOD;

  /*USER PARAMS:*/
  REAL nu;               /* GHz */
  REAL magnetic_field;   /* Gauss */
  REAL electron_density; /* 1/cc */
  REAL observer_angle;   /* rad */  
  int    distribution;     
  int    polarization; 
  int    mode;             /*Emissivity or Absorptivity*/
  REAL gamma_cutoff;
  /* Options for fits */
  int approximate;         /*Use approximate Bessel functions for speed*/
  int dexter_fit;

  /*Thermal distribution parameters*/
  REAL theta_e;

  /*power law parameters*/
  REAL power_law_p;
  REAL gamma_min;
  REAL gamma_max;

  /*kappa distribution parameters*/
  REAL kappa;
  REAL kappa_width;

  /*Choose if n-space peak is known, or if it must be found adaptively */
  int use_n_peak;
  REAL (*n_peak)(struct parameters *);

  /*Set distribution_function */
  REAL (*distribution_function)(REAL gamma, struct parameters *);

  /*analytic_differential_of_f, which can be used as a test of the 
    numerical differential_of_f */
  REAL (*analytic_differential)(REAL gamma, struct parameters *);

  int stokes_v_switch;

  char *error_message; /* if not NULL, records source of error in current calculation */

  /*susceptibility tensor paramsS */
  REAL gamma;
  REAL epsilon0;
  REAL epsilon;
  REAL omega;
  REAL omega_c;
  REAL omega_p;
  REAL real;
  REAL (*tau_integrand)(REAL, void * parameters);
  REAL (*gamma_integrand)(REAL, void * parameters);
};

struct parametersGSL
{
  struct parameters params;
  REAL n;
};

void setConstParams(struct parameters *params);
REAL get_nu_c(struct parameters params);
REAL get_omega_p(struct parameters params);

#endif /* SYMPHONY_PARAMS_H_ */
