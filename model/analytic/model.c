#include "model.h"

#include "decs.h"
#include "coordinates.h"
#include "geometry.h"
#include "radiation.h"
#include "par.h"

// Globals we're in charge of
REAL M_unit;
REAL L_unit;
REAL T_unit;
REAL U_unit;
REAL B_unit;
REAL Te_unit;
// This one is useful
REAL RHO_unit;

// Model parameters: public
REAL rmax_geo = 1000.0;
// Model parameters: private
static REAL rmin_geo = 0;
static REAL MBH_solar = 4.e6;
static REAL RHO_unit_in = 0;

static REAL MBH;
static int model;

// e.g. parameterization from GRRT paper
REAL A, alpha, height, l0, freqcgs;
REAL r_isco;

/**
 * This is a template for analytic problems, which consist of prescription:
 * X,K -> 4-vectors Ucon/cov, Bcon/cov
 * And either:
 * X,K -> emission coefficients jS, alphaS, rhoS
 * or:
 * X,K -> e- density and temperature ne, Thetae
 */

// Forward declarations for non-public functions
void set_units();

//// INITIALIZATION: Functions called from elsewhere in ipole ////

/**
 * This function is called for each word/value pair ipole encounters,
 * either from a parfile or the command line.
 * You can define new pairs here, skim what you need of ipole's defaults, or whatever
 * 
 * ipole will not warn on unspecified parameters. Have good defaults (set on declaration)
 */
void try_set_model_parameter(const char *word, const char *value)
{
  // Test the given word against our parameters' names,
  // and if it matches set the corresponding global
  set_by_word_val(word, value, "model", &model, TYPE_INT);
  set_by_word_val(word, value, "MBH", &MBH_solar, TYPE_DBL);
  set_by_word_val(word, value, "rho_unit", &RHO_unit_in, TYPE_DBL);
  // TODO NEED to move this into main parameters
  set_by_word_val(word, value, "rmax_geo", &rmax_geo, TYPE_DBL);
  set_by_word_val(word, value, "rmin_geo", &rmin_geo, TYPE_DBL);

  // Normal ipole pulls this, but we also need it for the GRRT problems
  // and this is easier than grabbing it from the 'params' struct
  set_by_word_val(word, value, "freqcgs", &freqcgs, TYPE_DBL);
}

/**
 * Initialization takes boundary times, for slow light.  Most analytic models won't use them.
 */
void init_model(REAL *tA, REAL *tB)
{
  // Set all the geometry globals we need
  // TODO do this in geometry?  Deal with model/geom interface...
  use_eKS_internal = 0;
  metric = METRIC_MKS;
  hslope = 1.0;

  if (model == 1) {
    A = 0;
    alpha = -3;
    height = 0;
    l0 = 0;
    a = 0.9;
  } else if (model == 2) {
    A = 0;
    alpha = -2;
    height = 0;
    l0 = 1;
    a = 0;
  } else if (model == 3) {
    A = 0;
    alpha = 0;
    height = 10./3;
    l0 = 1;
    a = 0.9;
  } else if (model == 4) {
    A = 1.e5;
    alpha = 0;
    height = 10./3;
    l0 = 1;
    a = 0.9;
  } else if (model == 5) {
    A = 1.e6;
    alpha = 0;
    height = 100./3;
    l0 = 1;
    a = 0.9;
  }

  // We already set stuff from parameters, so set_units here
  set_units();

  printf("Running analytic model %d:\nMBH: %Lg\na: %Lg\nRh: %Lg\nR_isco: %Lg\n\n", model, MBH, a, Rh, r_isco);
  printf("A: %Lg\nalpha: %Lg\nh: %Lg\nl0: %Lg\n\n", A, alpha, height, l0);
}

void set_units()
{
  // Derive units we actually need
  MBH = MBH_solar * MSUN;
  L_unit = GNEWT * MBH / (CL * CL);
  T_unit = L_unit / CL;
  // Could also set M_unit here from RHO_unit
  if (RHO_unit_in > 0) {
    RHO_unit = RHO_unit_in;
  } else {
    RHO_unit = 3.e-18;
  }
  B_unit = CL * sqrt(4.*M_PI*RHO_unit);

  // Set all the geometry for coordinates.c
  // TODO function like initialize_coordinates, that makes sure these are all set.
  R0 = 0.;
  Rh = 1 + sqrt(1. - a * a);
  Rin = Rh;
  Rout = 1000.0;
  // Limit rmax_geo?

  REAL z1 = 1. + pow(1. - a * a, 1. / 3.) * (pow(1. + a, 1. / 3.) + pow(1. - a, 1. / 3.));
  REAL z2 = sqrt(3. * a * a + z1 * z1);
  r_isco = 3. + z2 - copysign(sqrt((3. - z1) * (3. + z1 + 2. * z2)), a);
  startx[0] = 0.0;
  startx[1] = log(Rin);
  startx[2] = 0.0;
  startx[3] = 0.0;
  stopx[0] = 0.0;
  stopx[1] = log(Rout);
  stopx[2] = 1.0;
  stopx[3] = 2*M_PI;
}

void output_hdf5()
{
  hdf5_set_directory("/header/");
  double zero = 0;
  hdf5_write_single_val(&zero, "t", H5T_DOUBLE);
  hdf5_write_single_val(&a, "a", H5T_DOUBLE);

  hdf5_write_single_val(&model, "model", H5T_STD_I32LE);
  hdf5_write_single_val(&A, "A", H5T_DOUBLE);
  hdf5_write_single_val(&alpha, "alpha", H5T_DOUBLE);
  hdf5_write_single_val(&height, "height", H5T_DOUBLE);
  hdf5_write_single_val(&l0, "l0", H5T_DOUBLE);

  hdf5_make_directory("units");
  hdf5_set_directory("/header/units/");
  hdf5_write_single_val(&zero, "M_unit", H5T_DOUBLE);
  hdf5_write_single_val(&L_unit, "L_unit", H5T_DOUBLE);
  hdf5_write_single_val(&T_unit, "T_unit", H5T_DOUBLE);

  hdf5_set_directory("/");
}

//// INTERFACE: Functions called from elsewhere in ipole ////
REAL get_model_ne(REAL X[NDIM])
{
  // Matter model defined in Gold et al 2020 section 3
  REAL r, th;
  bl_coord(X, &r, &th);
  REAL n_exp = 1./2 * (pow(r/10, 2) + pow(height * cos(th), 2));
  // Cutoff when result will be ~0
  return ( n_exp < 200 ) ? RHO_unit * exp(-n_exp) : 0;
}

void get_model_jk(REAL X[NDIM], REAL Kcon[NDIM], REAL *jnuinv, REAL *knuinv)
{
  // Emission model defined in Gold et al 2020 section 3
  REAL n = get_model_ne(X);

  REAL Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
  get_model_fourv(X, Kcon, Ucon, Ucov, Bcon, Bcov);
  REAL nu = get_fluid_nu(Kcon, Ucov);

  *jnuinv = fmax( n * pow(nu / freqcgs, -alpha) / pow(nu, 2), 0);
  *knuinv = fmax( (A * n * pow(nu / freqcgs, -(2.5 + alpha)) + 1.e-54) * nu, 0);
}

void get_model_jar(REAL X[NDIM], REAL Kcon[NDIM],
    REAL *jI, REAL *jQ, REAL *jU, REAL *jV,
    REAL *aI, REAL *aQ, REAL *aU, REAL *aV,
    REAL *rQ, REAL *rU, REAL *rV)
{
  // Define a model here relating X,K -> j_S, alpha_S, rho_S
  // (and below relating X,K -> u,B 4-vectors)
  // ipole will do the rest

  // Just take the unpolarized emissivity and absorptivity as I,
  // and set the rest to zero
  // Of course, you can be more elaborate
  REAL j, k;
  get_model_jk(X, Kcon, &j, &k);

  *jI = j;
  *jQ = 0;
  *jU = 0;
  *jV = 0;

  *aI = k;
  *aQ = 0;
  *aU = 0;
  *aV = 0;

  *rQ = 0;
  *rU = 0;
  *rV = 0;

  return;
}

void get_model_fourv(REAL X[NDIM], REAL Kcon[NDIM], REAL Ucon[NDIM], REAL Ucov[NDIM],
                     REAL Bcon[NDIM], REAL Bcov[NDIM])
{
  REAL r, th;
  bl_coord(X, &r, &th);
  // Note these quantities from Gold et al are in BL!
  // We could have converted the problem to KS, but instead we did this
  REAL R = r * sin(th);
  REAL l = (l0 / (1 + R)) * pow(R, 1 + 0.5);

  // Metrics: BL
  REAL bl_gcov[NDIM][NDIM], bl_gcon[NDIM][NDIM];
  gcov_bl(r, th, bl_gcov);
  gcon_func(bl_gcov, bl_gcon);
  // Native
  REAL gcov[NDIM][NDIM], gcon[NDIM][NDIM];
  gcov_func(X, gcov);
  gcon_func(gcov, gcon);

  // Get the normal observer velocity for Ucon/Ucov, in BL coordinates
  REAL bl_Ucov[NDIM];
  REAL ubar = sqrt(-1. / (bl_gcon[0][0] - 2. * bl_gcon[0][3] * l
                  + bl_gcon[3][3] * l * l));
  bl_Ucov[0] = -ubar;
  bl_Ucov[1] = 0.;
  bl_Ucov[2] = 0.;
  bl_Ucov[3] = l * ubar;

  REAL bl_Ucon[NDIM];
  flip_index(bl_Ucov, bl_gcon, bl_Ucon);

  // Transform to KS coordinates,
  REAL ks_Ucon[NDIM];
  bl_to_ks(X, bl_Ucon, ks_Ucon);
  // then to our coordinates,
  vec_from_ks(X, ks_Ucon, Ucon);

  // and grab Ucov
  flip_index(Ucon, gcov, Ucov);

  // ...or don't do any of that
  // REAL ubar = sqrt(-1. / (gcon[0][0] - 2. * gcon[0][3] * l
  //                 + gcon[3][3] * l * l));
  // Ucov[0] = -ubar;
  // Ucov[1] = 0.;
  // Ucov[2] = 0.;
  // Ucov[3] = l * ubar;
  // flip_index(Ucov, gcon, Ucon);

  // This model defines no field in emission, but the field is used for making
  // tetrads so we want it consistent
  Bcon[0] = 0;
  Bcon[1] = 0;
  Bcon[2] = 1;
  Bcon[3] = 0;
  flip_index(Bcon, gcov, Bcov);
}

/**
 * This problem defines no field in emission, but we want to control ipole's worst
 * tendencies when making tetrads.  This will return a correct value even for a
 * possible fluid/field model later, too.
 */
REAL get_model_b(REAL X[NDIM])
{
  REAL Ucon[NDIM],Bcon[NDIM];
  REAL Ucov[NDIM],Bcov[NDIM];
  REAL Kcon[NDIM] = {0}; // TODO interface change if we ever need a real one here
  get_model_fourv(X, Kcon, Ucon, Ucov, Bcon, Bcov);
  return sqrt(Bcon[0]*Bcov[0] + Bcon[1]*Bcov[1] + Bcon[2]*Bcov[2] + Bcon[3]*Bcov[3]) * B_unit;
}

int radiating_region(REAL X[NDIM])
{
  // If you don't want conditionals in get_model_jar, 
  // you can control here where the coefficients are applied
  REAL r, th;
  bl_coord(X, &r, &th);
  return r > Rh + 0.0001 && r > rmin_geo && r < 1000.0;
}

//// STUBS: Functions for normal models which we don't use ////
// Define these to specify a fluid model: e- density/temperature for
// synchrotron radiation based on an energy distribution
REAL get_model_thetae(REAL X[NDIM]) {return 0;}
void get_model_powerlaw_vals(REAL X[NDIM], REAL *p, REAL *n,
          REAL *gamma_min, REAL *gamma_max, REAL *gamma_cut) {return;}

// This is only called for trace file output, and doesn't really apply to analytic models
void get_model_primitives(REAL X[NDIM], REAL *p) {return;}
