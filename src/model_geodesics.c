
/*
 * model-dependent routines for integrating geodesics,
 * including:
 *
 * stop_backward_integration
 * stepsize
 */

#include "model_geodesics.h"

#include "coordinates.h"
#include "decs.h"
#include "geodesics.h"
#include "geometry.h"
#include "model.h"
#include "model_tetrads.h"
#include "tetrads.h"

#include "debug_tools.h"

// TODO pick one or runtime
#define STEP_STRICT_MIN 0

/**
 * Trace a single geodesic
 * Takes a starting location and frequency in *units of electron mass frequency*,
 * as well as an array of trajectory location structs of_traj, size MAXNSTEP
 * 
 * TODO this routine definitely doesn't initialize traj[0].
 * 
 * Note convention changes for integrating backwards:
 * * Xhalf & Kconhalf trail X & Kcon in this function, and thereafter lead them
 * * In this function, dl is the length of the step *to* point N;
 *   afterward it is *from* point N onward
 */
int trace_geodesic(REAL Xi[NDIM], REAL Kconi[NDIM], struct of_traj *traj, REAL eps, int step_max)
{
  //fprintf(stderr, "Begin trace geodesic");
  REAL X[NDIM], Kcon[NDIM];
  REAL Xhalf[NDIM], Kconhalf[NDIM];

  // Initialize the running values and save the first step
  // Note "half" values *trail* when integrating away from camera, so they don't
  // make sense to record -- we just initialize them for safety
  MULOOP {
    X[mu] = Xi[mu];
    Xhalf[mu] = Xi[mu];
    Kcon[mu] = Kconi[mu];
    Kconhalf[mu] = Kconi[mu];

    traj[0].X[mu] = Xi[mu];
    traj[0].Kcon[mu] = Kconi[mu];
  }

  int nstep = 0;

  // Integrate backwards
  while ( (!stop_backward_integration(X, Xhalf, Kcon)) && (nstep < step_max - 1) ) {
    nstep++;

    /* This stepsize function can be troublesome inside of R = 2M,
       and should be used cautiously in this region. */
    REAL dl = stepsize(X, Kcon, eps);

    /* Geodesics in ipole are integrated using
     * dx^\mu/d\lambda = k^\mu
     * The positions x^mu are in simulation units, since different
     * coordinates sometimes have different units (e.g. x^r, x^\theta)
     *
     * The convention we have adopted is that:
     * E = -u^\mu k_\mu,
     * which is always photon energy measured by an observer with four-velocity u^\mu,
     * is in units of *electron rest-mass energy*.
     * This implies that dl is *not* in cgs units, but in weird hybrid units.
     * This line sets dl to be in cgs units.
     */

#if INTEGRATOR_TEST
    traj[nstep].dl = dl;
#else
    traj[nstep].dl = dl * L_unit * HPL / (ME * CL * CL);
#endif

    // To print each point (TODO option?)
    // print_vector("X", X);
    // print_vector("Kcon", Kcon);
    // fprintf(stderr, "dl: %Lf\n", dl);
    // REAL r, th;
    // bl_coord(X, &r, &th);
    // fprintf(stderr, "BL r, th: %Lg, %Lg\n", r, th);
    // REAL Gcov[NDIM][NDIM];
    // gcov_func(X, Gcov);
    // print_matrix("Gcov", Gcov);
    // getchar();

    /* move photon one step backwards, the procecure updates X
       and Kcon full step and returns also values in the middle */
    push_photon(X, Kcon, -dl, Xhalf, Kconhalf);

    // Set the new position and wavevector
    MULOOP {
      traj[nstep].X[mu] = X[mu];
      traj[nstep].Kcon[mu] = Kcon[mu];
      traj[nstep].Xhalf[mu] = Xhalf[mu];
      traj[nstep].Kconhalf[mu] = Kconhalf[mu];
    }
  }

  //fprintf(stderr, "End trace geodesic");
  return nstep;
}

/*
 * Initialize a geodesic from the camera
 * This takes the parameters struct directly since most of them are
 * camera parameters anyway
 */
void init_XK(long int i, long int j, int nx, int ny, REAL Xcam[NDIM],
             Params params, REAL fovx, REAL fovy,
             REAL X[NDIM], REAL Kcon[NDIM])
{
  REAL Econ[NDIM][NDIM];
  REAL Ecov[NDIM][NDIM];
  REAL Kcon_tetrad[NDIM];

  if (params.old_centering) {
    make_camera_tetrad_old(Xcam, Econ, Ecov);
  } else {
    make_camera_tetrad(Xcam, Econ, Ecov);
  }

  // Construct outgoing wavevectors
  // xoff: allow arbitrary offset for e.g. ML training imgs
  // +0.5: project geodesics from px centers
  // xoff/yoff are separated to keep consistent behavior between refinement levels
  REAL dxoff = (i+0.5+params.xoff-0.01)/params.nx - 0.5;
  REAL dyoff = (j+0.5+params.yoff)/params.ny - 0.5;
  Kcon_tetrad[0] = 0.;
  Kcon_tetrad[1] = (dxoff*cos(params.rotcam) - dyoff*sin(params.rotcam)) * fovx;
  Kcon_tetrad[2] = (dxoff*sin(params.rotcam) + dyoff*cos(params.rotcam)) * fovy;
  Kcon_tetrad[3] = 1.;

  /* normalize */
  null_normalize(Kcon_tetrad, 1.);

  /* translate into coordinate frame */
  tetrad_to_coordinate(Econ, Kcon_tetrad, Kcon);

  /* set position */
  // TODO this doesn't initialize X0 properly
  MULOOP X[mu] = Xcam[mu];
}


/* condition for stopping the backwards-in-lambda
   integration of the photon geodesic */

int stop_backward_integration(REAL X[NDIM], REAL Xhalf[NDIM], REAL Kcon[NDIM])
{
  // The opaque thin disk adds a stop condidion: we don't bother integrating more than 2 steps
  // beyond the midplane
#if THIN_DISK
  static int n_left = -1;
#pragma omp threadprivate(n_left)
#endif

  // Necessary geometric stop conditions
  REAL r, th;
  bl_coord(X, &r, &th);
  if ((r > rmax_geo && Kcon[1] < 0.) || // Stop either beyond rmax_geo
      r < (Rh + 0.0001)) { // Or right near the horizon
#if THIN_DISK
    // If we stopped during the thin disk timer, remember to reset it!
    n_left = -1;
#endif
    return (1);
  }

  // Additional stop condition for thin disks: the disk is opaque
#if THIN_DISK
  if (thindisk_region(X, Xhalf) && n_left < 0) { // Set timer when we reach disk
    n_left = 2;
    return 0;
  } else if (n_left > 0) { // Or decrement the timer if it's set
    n_left--;
    return 0;
  } else if (n_left == 0) { // If timer is 0, stop and reset
    n_left = -1;
    return 1;
  }
#endif

  return (0);
}


/**
 * Choose stepsize according to inverse Kcon, dramatically decreasing the step
 * toward the coordinate pole and EH.
 * 
 * Use the sum of inverses by default; the strict minimum seems to occasionally
 * overstep even for small eps
 * 
 * TODO this is the geometry step but dictates physics as well.
 * Optionally skip geometry steps near the pole for accuracy
 */
REAL stepsize(REAL X[NDIM], REAL Kcon[NDIM], REAL eps)
{
  REAL dl;
  REAL deh = fmin(fabs(X[1] - startx[1]), 0.1);
  REAL dlx1 = eps * (10*deh) / (fabs(Kcon[1]) + SMALL*SMALL);

  // Make the step cautious near the pole, improving accuracy of Stokes U
  REAL cut = 0.02;
  REAL lx2 = stopx[2] - startx[2];
  REAL dpole = fmin(fabs(X[2] / lx2), fabs((stopx[2] - X[2]) / lx2));
  REAL d2fac = (dpole < cut) ? dpole/3 : fmin(cut/3 + (dpole-cut)*10., 1);
  REAL dlx2 = eps * d2fac / (fabs(Kcon[2]) + SMALL*SMALL);

  REAL dlx3 = eps / (fabs(Kcon[3]) + SMALL*SMALL);

  if (STEP_STRICT_MIN) {
    dl = fmin(fmin(dlx1, dlx2), dlx3);
  } else {
    REAL idlx1 = 1./(fabs(dlx1) + SMALL*SMALL) ;
    REAL idlx2 = 1./(fabs(dlx2) + SMALL*SMALL) ;
    REAL idlx3 = 1./(fabs(dlx3) + SMALL*SMALL) ;

    dl = 1. / (idlx1 + idlx2 + idlx3) ;
  }

  return dl;
}

