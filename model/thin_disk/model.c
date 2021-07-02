// Standard NT73 thin disk model
// Pretty much cribbed from grtrans

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
REAL RHO_unit;
REAL U_unit;
REAL B_unit;
REAL Te_unit;

// Model parameters: public
REAL rmax_geo;
// Model parameters: private
static REAL Mdot = 0.01;
static REAL MBH_solar = 10.;
static REAL MBH;

// Thin-disk tables
REAL *ch_mu, *ch_I, *ch_delta;
// Thin-disk constants
static REAL T0, Mdotedd, r_isco;

// Forward declarations for non-public functions
void set_units();
void thindisk_vals(REAL r, REAL *T, REAL *omega);
void calc_polvec(REAL X[NDIM], REAL Kcon[NDIM], REAL a, REAL fourf[NDIM]);
REAL krolikc(REAL r, REAL a);
// Emission
void fbbpolemis(REAL nu, REAL T, REAL cosne, REAL *SI, REAL *SQ);
REAL bnu(REAL nu, REAL T);
// Chandrasekhar table stuff
void load_chandra_tab24();
void interp_chandra(REAL mu, REAL *i, REAL *del);


void try_set_model_parameter(const char *word, const char *value)
{
  // Test the given word against our parameters' names,
  // and if it matches set the corresponding global
  set_by_word_val(word, value, "MBH", &MBH_solar, TYPE_DBL);

  set_by_word_val(word, value, "a", &a, TYPE_DBL);
  set_by_word_val(word, value, "Mdot", &Mdot, TYPE_DBL);
}

void init_model(REAL *tA, REAL *tB)
{
  load_chandra_tab24();

  // Set all the geometry globals we need
  // TODO do this in geometry?  Deal with model/geom interface...
  use_eKS_internal = 1;
  metric = 0; // Doesn't matter due to above
  hslope = 1.0;
  // Needed for camera rootfinding
  // TODO standard geometry init that handles this...
  startx[2] = 0.;
  stopx[2] = 1;

  // We already set stuff from parameters, so set_units here
  set_units();

  printf("Running thin disk:\nMBH: %Lg\nMdot: %Lg\na: %Lg\n\n", MBH, Mdot, a);
}

void set_units()
{
  // Convert parameters to consistent CGS units
  MBH = MBH_solar * MSUN;
  // Note definition of Mdotedd w/ 10% efficiency
  Mdotedd = 4. * M_PI * GNEWT * MBH * MP / 0.1 / CL / SIGMA_THOMSON;
  Mdot *= Mdotedd;

  // Derive everything else
  M_unit = Mdot;
  L_unit = GNEWT * MBH / (CL * CL);
  T_unit = L_unit / CL;
  RHO_unit = M_unit / pow(L_unit, 3.);
  U_unit = RHO_unit * CL * CL;
  B_unit = CL * sqrt(4. * M_PI * RHO_unit);

  // Set all the geometry
  // TODO function like initialize_coordinates, that makes sure these are all set.
  R0 = 0.;
  Rh = 1 + sqrt(1. - a * a);
  REAL z1 = 1. + pow(1. - a * a, 1. / 3.) * (pow(1. + a, 1. / 3.) + pow(1. - a, 1. / 3.));
  REAL z2 = sqrt(3. * a * a + z1 * z1);
  r_isco = 3. + z2 - copysign(sqrt((3. - z1) * (3. + z1 + 2. * z2)), a);
  Rin = Rh;
  Rout = 100.0;
  rmax_geo = fmin(1000., Rout);
  startx[0] = 0.0;
  startx[1] = log(Rin);
  startx[2] = 0.0;
  startx[3] = 0.0;
  stopx[0] = 0.0;
  stopx[1] = log(Rout);
  stopx[2] = 1.0;
  stopx[3] = 2*M_PI;

  // Some precomputation.  TODO should this be Mdot?
  T0 = pow(3.0 / 8.0 / M_PI * GNEWT * MBH * M_unit / pow(L_unit, 3) / SIG, 1. / 4.);

  fprintf(stderr, "L,T,M units: %Lg [cm] %Lg [s] %Lg [g]\n", L_unit, T_unit, M_unit);
  fprintf(stderr, "rho,u,B units: %Lg [g cm^-3] %Lg [g cm^-1 s^-2] %Lg [G] \n", RHO_unit, U_unit, B_unit);
  fprintf(stderr, "Rh, Rin, Rout, r_isco, T0: %Lg %Lg %Lg %Lg %Lg\n", Rh, Rin, Rout, r_isco, T0);
}

void output_hdf5()
{
  hdf5_set_directory("/header/");
  REAL zero = 0;
  hdf5_write_single_val(&zero, "t", H5T_DOUBLE);
  hdf5_write_single_val(&a, "a", H5T_DOUBLE);
  hdf5_write_single_val(&Mdot, "Mdot", H5T_DOUBLE);

  hdf5_make_directory("units");
  hdf5_set_directory("/header/units/");
  hdf5_write_single_val(&L_unit, "L_unit", H5T_DOUBLE);
  hdf5_write_single_val(&M_unit, "M_unit", H5T_DOUBLE);
  hdf5_write_single_val(&T_unit, "T_unit", H5T_DOUBLE);
  //hdf5_write_single_val(&Te_unit, "Thetae_unit", H5T_DOUBLE);

  hdf5_set_directory("/");
}

//// PUBLIC INTERFACE ////
void get_model_stokes(REAL X[NDIM], REAL Kcon[NDIM], REAL *SI,
                      REAL *SQ, REAL *SU, REAL *SV)
{
  REAL r, th;
  bl_coord(X, &r, &th);

  if (r > Rh) {
    REAL T, omega;
    thindisk_vals(r, &T, &omega);

    REAL Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
    get_model_fourv(X, Kcon, Ucon, Ucov, Bcon, Bcov);

    // Recall "B" was set by calc_polvec
    REAL mu = fabs(cos(get_bk_angle(X, Kcon, Ucov, Bcon, Bcov)));
    REAL nu = get_fluid_nu(Kcon, Ucov);

    fbbpolemis(nu, T, mu, SI, SQ);

  } else {
    *SI = 0;
    *SQ = 0;
  }

  *SU = 0;
  *SV = 0;
}

void get_model_i(REAL X[NDIM], REAL Kcon[NDIM], REAL *SI)
{
  // Unpolarized emission is calculated no differently for this model
  REAL DQ, DU, DV;
  get_model_stokes(X, Kcon, SI, &DQ, &DU, &DV);
}

int thindisk_region(REAL Xi[NDIM], REAL Xf[NDIM])
{
  REAL ri, thi, rf, thf;
  bl_coord(Xi, &ri, &thi);
  bl_coord(Xf, &rf, &thf);
  // Set the intensity whenever and exactly when we cross the disk, outside horizon
  int midplane = (sign(thi - M_PI_2) != sign(thf - M_PI_2));
  int em_region = rf > r_isco && rf < Rout;
  return midplane && em_region;

  // Less restrictive
  //return (fabs(fabs(thf) - M_PI_2) < 0.01);
}

void get_model_fourv(REAL X[NDIM], REAL Kcon[NDIM], REAL Ucon[NDIM],
                       REAL Ucov[NDIM], REAL Bcon[NDIM], REAL Bcov[NDIM])
{
  REAL r, th, T, omega;
  REAL gcov[NDIM][NDIM];

  // Get native metric
  gcov_func(X, gcov);

  // Then get some stuff in BL for Ucon
  bl_coord(X, &r, &th);

  thindisk_vals(r, &T, &omega);

  // normal observer velocity for Ucon/Ucov
  Ucon[0] =
      sqrt(-1. / (gcov[0][0] + 2. * gcov[0][3] * omega
                  + gcov[3][3] * omega * omega));
  Ucon[1] = 0.;
  Ucon[2] = 0.;
  Ucon[3] = omega * Ucon[0];

  flip_index(Ucon, gcov, Ucov);

  // B is handled in native coordinates here
  calc_polvec(X, Kcon, a, Bcon);

  // Flip B
  flip_index(Bcon, gcov, Bcov);
}

//// SUPPORT: Thin Disk functions ////

// Only supports midplane!
void thindisk_vals(REAL r, REAL *T, REAL *omega)
{
  REAL b, kc, d, ar, lc, hc, om;

  REAL th = M_PI_2;

  b = 1. - 3. / r + 2. * a / pow(r, 3. / 2.);
  kc = krolikc(r, a);
  d = r * r - 2. * r + a * a;
  lc = (r_isco * r_isco - 2. * a * sqrt(r_isco) + a * a)
      / (pow(r_isco, 1.5) - 2. * sqrt(r_isco) + a);
  hc = (2. * r - a * lc) / d;

  ar = pow(r * r + a * a, 2.) - a * a * d * pow(sin(th), 2.);
  om = 2. * a * r / ar;

  // Start the disk at r_isco, the marginally stable orbit which N-K take as an inner boundary condition.
  // End it eventually.
  if (r > r_isco) {
    *omega = fmax(1. / (pow(r, 3. / 2.) + a), om);
  } else {
    *omega = fmax((lc + a * hc) / (r * r + 2. * r * (1. + hc)), om);
  }

  if (r > r_isco && r < Rout) {
    *T = T0 * pow(kc / b / pow(r, 3), 1. / 4.);
  } else {
    *T = T0 / 1e5;
  }
}

REAL krolikc(REAL r, REAL a)
{
  REAL y = sqrt(r);
  REAL yms = sqrt(r_isco);
  REAL y1 = 2. * cos(1. / 3. * (acos(a) - M_PI));
  REAL y2 = 2. * cos(1. / 3. * (acos(a) + M_PI));
  REAL y3 = -2. * cos(1. / 3. * acos(a));
  REAL arg1 = 3. * a / (2. * y);
  REAL arg2 = 3. * pow(y1 - a, 2) / (y * y1 * (y1 - y2) * (y1 - y3));
  REAL arg3 = 3. * pow(y2 - a, 2) / (y * y2 * (y2 - y1) * (y2 - y3));
  REAL arg4 = 3. * pow(y3 - a, 2) / (y * y3 * (y3 - y1) * (y3 - y2));

  return 1. - yms / y - arg1 * log(y / yms) - arg2 * log((y - y1) / (yms - y1))
      - arg3 * log((y - y2) / (yms - y2)) - arg4 * log((y - y3) / (yms - y3));
}

//// SUPPORT: Emissivities ////

/*
 * Set photon wavevector for each radial zone of thin disk (polarized emission)
 */
void fbbpolemis(REAL nu, REAL T, REAL cosne, REAL *SI, REAL *SQ)
{
  REAL f = 1.8;
  *SI = pow(f, -4.) * bnu(nu, T * f);

  // assumes Chandrasekhar electron scattering from semi-infinite atmosphere
  REAL interpI, interpdel;
  interp_chandra(cosne, &interpI, &interpdel);
  *SI = *SI * interpI;
  *SQ = *SI * interpdel;

  // Return invariant intensity & polarization
  *SI = *SI / (nu*nu*nu);
  *SQ = *SQ / (nu*nu*nu);
}

// Blackbody function B_nu(theta_e)
// c.f. Bnu_inv, same function but different interface
REAL bnu(REAL nu, REAL T) {
  return 2 * HPL * nu*nu*nu / (CL * CL) / (exp(HPL * nu / (KBOL * T)) - 1);
}

//// SUPPORT: Tetrads ////

// This sure is a vector.
void calc_polvec(REAL X[NDIM], REAL Kcon[NDIM], REAL a, REAL fourf[NDIM])
{
  REAL fourf_bl[NDIM], fourf_ks[NDIM];
  fourf_bl[0] = 0;
  fourf_bl[1] = 0;
  fourf_bl[2] = 1;
  fourf_bl[3] = 0;

  // Then transform to KS and to eKS
  bl_to_ks(X, fourf_bl, fourf_ks);
  vec_to_ks(X, fourf_ks, fourf);

  // Now normalize
  REAL gcov[NDIM][NDIM], fourf_cov[NDIM];
  gcov_func(X, gcov);
  flip_index(fourf, gcov, fourf_cov);
  REAL normf = sqrt(fourf[0] * fourf_cov[0] + fourf[1] * fourf_cov[1] + fourf[2] * fourf_cov[2]
      + fourf[3] * fourf_cov[3]);

  MULOOP fourf[mu] /= normf;
}

// Functions for dealing with Chandrasekhar 1960 table 24,
// i.e. polarized emission via scattering from a semi-infinite atmosphere
void load_chandra_tab24()
{
  char *fname = "ch24_vals.txt";
  FILE *vals;
  vals = fopen(fname, "r");
  if (vals == 0){
    fprintf(stderr, "Error reading file %s!\n\n", fname);
    exit(-1);
  }

  int npts = 21; // TODO read at runtime?  Set up for it already...
  ch_mu = calloc(sizeof(REAL), npts);
  ch_I = calloc(sizeof(REAL), npts);
  ch_delta = calloc(sizeof(REAL), npts);
  int gcc_stahp;
  for (int i = 0; i < npts; ++i) {
    gcc_stahp = fscanf(vals, "%lf", &ch_mu[i]);
    gcc_stahp = fscanf(vals, "%lf", &ch_I[i]);
    gcc_stahp = fscanf(vals, "%lf", &ch_delta[i]);
    //fprintf(stderr, "TABLE: %Lf %Lf %Lf\n", ch_mu[i], ch_I[i], ch_delta[i]);
  }
  fclose(vals);
}

REAL get_weight(REAL *xx, REAL x, int *jlo)
{
  //Get the value _before_ x in the table
  while (xx[*jlo] < x) {
    ++*jlo;
  }
  --*jlo;

  // Return weight for table values jlo, jlo+1
  return (x - xx[*jlo]) / (xx[*jlo + 1] - xx[*jlo]);
}

void interp_chandra(REAL mu, REAL *i, REAL *del)
{
  int indx = 0;
  REAL weight = get_weight(ch_mu, mu, &indx);
  *i = (1. - weight) * ch_I[indx] + weight * ch_I[indx + 1];
  *del = (1. - weight) * ch_delta[indx] + weight * ch_delta[indx + 1];
}


//// STUBS: Functions for normal models which we don't use ////
int radiating_region(REAL X[NDIM]) {return 0;}
REAL get_model_thetae(REAL X[NDIM]) {return 0;}
REAL get_model_b(REAL X[NDIM]) {return 0;}
REAL get_model_ne(REAL X[NDIM]) {return 0;}
void get_model_primitives(REAL X[NDIM], REAL *p) {return;}
void get_model_powerlaw_vals(REAL X[NDIM], REAL *p, REAL *n,
          REAL *gamma_min, REAL *gamma_max, REAL *gamma_cut) {return;}
// In case we want to mess with emissivities directly
void get_model_jar(REAL X[NDIM], REAL Kcon[NDIM],
    REAL *jI, REAL *jQ, REAL *jU, REAL *jV,
    REAL *aI, REAL *aQ, REAL *aU, REAL *aV,
    REAL *rQ, REAL *rU, REAL *rV) {return;}
void get_model_jk(REAL X[NDIM], REAL Kcon[NDIM], REAL *jnuinv, REAL *knuinv) {return;}