#include "decs.h"

#include "model.h"
#include "model_geodesics.h"
#include "model_radiation.h"
#include "model_tetrads.h"

#include "radiation.h"
#include "coordinates.h"
#include "debug_tools.h"
#include "tetrads.h"
#include "geometry.h"
#include "geodesics.h"
#include "image.h"
#include "io.h"
#include "ipolarray.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <omp.h>

// Print a backtrace of a certain pixel, useful for debugging
// Set to -1 to disable
#define DIAG_PX_I -1
#define DIAG_PX_J -1

// Some useful blocks of code to re-use
// Note the difference between "int nx,ny" in get_pixel and "long int nx,ny" in save_pixel
// explained in parameter parsing.
void get_pixel(size_t i, size_t j, int nx, int ny, REAL Xcam[NDIM], Params params,
               REAL fovx, REAL fovy, REAL freq, int only_intensity, REAL scale,
               REAL *Intensity, REAL *Is, REAL *Qs, REAL *Us, REAL *Vs,
               REAL *Tau, REAL *tauF);
void save_pixel(REAL *image, REAL *imageS, REAL *taus, size_t i, size_t j, size_t nx, size_t ny, int only_unpol,
                REAL Intensity, REAL Is, REAL Qs, REAL Us, REAL Vs,
                REAL freqcgs, REAL Tau, REAL tauF);
void print_image_stats(REAL *image, REAL *imageS, size_t nx, size_t ny, Params params, REAL scale);
void save_pixelTransfer(REAL *image, REAL *imageS, REAL *taus,
                        size_t iold, size_t jold, size_t inew, size_t jnew, size_t nx, size_t ny, int only_intensity); //nearest neighbor saving
void lininterp4(REAL *image, REAL *imageS, REAL *taus, size_t i1, size_t j1,
                size_t i2,size_t j2, size_t i3, size_t j3, size_t i4, size_t j4, size_t inew, size_t jnew,
                size_t nx, size_t ny, int only_intensity); //linear interpolation for floater case 
void lininterp2(REAL *image, REAL *imageS, REAL *taus, size_t i1, size_t j1,
                size_t i2,size_t j2,size_t inew, size_t jnew, size_t nx, size_t ny, int only_intensity); //linear interpolation for same row or column case


// TODO use this more.  Are long=only or upscaled ops faster?
static inline size_t imgindex(size_t n, size_t i, size_t j, size_t nx, size_t ny) {return (n*nx + i)*ny + j;}

// global variables. TODO scope into main
static REAL tf = 0.;

Params params = { 0 };

int main(int argc, char *argv[]) 
{
  // motd
  fprintf(stderr, "%s. githash: %s\n", VERSION_STRING, xstr(VERSION));
  fprintf(stderr, "notes: %s\n\n", xstr(NOTES));

  // initialization
  REAL time = omp_get_wtime();

  REAL tA, tB; // for slow light
  REAL Xcam[NDIM];
  REAL freq, scale;
  REAL DX, DY, fovx, fovy;

#pragma omp parallel
  if (omp_get_thread_num() == 0) {
    fprintf(stderr, "nthreads = %d\n", omp_get_num_threads());
  }

  // load values from parameter file. handle all actual
  // model parameter comprehension in the model/* files
  load_par_from_argv(argc, argv, &params);

  // now that we've loaded all parameters, tell our model about
  // them and use init_model to load the first dump
  init_model(&tA, &tB);

  // Adaptive resolution option
  // nx, ny are the resolution at maximum refinement level
  // nx_min, ny_min are at coarsest level
  // TODO check for obvious BS here
  if (params.nx_min < 0) {
    params.nx_min = params.nx;
    params.ny_min = params.ny;
  }

  int refine_level = log2((params.nx - params.nx % 2)/(params.nx_min - params.nx_min % 2))+1;
  // INTERNAL SIZE.  If nx or ny is even, compute an extra row/column to use the 2^N+1 scheme
  size_t nx, ny, nxmin, nymin;
  if (refine_level > 1 && params.nx % 2 == 0) {
    nx = params.nx + 1;
    nxmin = params.nx_min + 1;
  } else {
    nx = params.nx;
    nxmin = params.nx_min;
  }
  if (refine_level > 1 && params.ny % 2 == 0) {
    ny = params.ny + 1;
    nymin = params.ny_min + 1;
  } else {
    ny = params.ny;
    nymin = params.ny_min;
  }

  // normalize frequency to electron rest-mass energy
  REAL freqcgs = params.freqcgs;
  freq = params.freqcgs * HPL / (ME * CL * CL);

  // Initialize the camera
  params.rotcam *= M_PI/180.;

  // translate to geodesic coordinates
  native_coord(params.rcam, params.thetacam, params.phicam, Xcam);
  fprintf(stderr, "Xcam[] = %Le %Le %Le %Le\n", Xcam[0], Xcam[1], Xcam[2], Xcam[3]);

  params.dsource *= PC;
  REAL Dsource = params.dsource; // Shorthand

  // set DX/DY using fov_dsource if possible, otherwise DX, otherwise old default
  REAL fov_to_d = Dsource / L_unit / MUAS_PER_RAD;

  if (params.fovx_dsource != 0.0) { // FOV was specified
    // Uncomment to be even more option lenient
    //if (params.fovy_dsource == 0.0) params.fovy_dsource = params.fovx_dsource;
  } else if (params.dx != 0.0) {
    //if (params.dy == 0.0) params.dy = params.dx;
    params.fovx_dsource = params.dx / fov_to_d;
    params.fovy_dsource = params.dy / fov_to_d;
  } else {
    fprintf(stderr, "No FOV was specified. Using default 160muas!\n");
    params.fovx_dsource = 160.;
    params.fovy_dsource = 160.;
  }
  DX = params.fovx_dsource * fov_to_d;
  DY = params.fovy_dsource * fov_to_d;
  params.dx = DX;
  params.dy = DY;

  // Set the *camera* fov values
  // We don't set these like other parameters, but output them for historical reasons
  fovx = DX / params.rcam;
  fovy = DY / params.rcam;

  scale = (DX * L_unit / params.nx) * (DY * L_unit / params.ny) / (Dsource * Dsource) / JY;
  fprintf(stderr,"intensity [cgs] to flux per pixel [Jy] conversion: %Lg\n",scale);
  fprintf(stderr,"Dsource: %Lg [cm]\n",Dsource);
  fprintf(stderr,"Dsource: %Lg [kpc]\n",Dsource/(1.e3*PC));
  fprintf(stderr,"FOVx, FOVy: %Lg %Lg [GM/c^2]\n",DX,DY);
  fprintf(stderr,"FOVx, FOVy: %Lg %Lg [rad]\n",DX*L_unit/Dsource,DY*L_unit/Dsource);
  fprintf(stderr,"FOVx, FOVy: %Lg %Lg [muas]\n",DX*L_unit/Dsource * MUAS_PER_RAD ,DY*L_unit/Dsource * MUAS_PER_RAD);
  if (refine_level > 1) {
    fprintf(stderr,"Resolution: %dx%d, refined up to %dx%d (%d levels)\n",
        params.nx_min, params.ny_min, params.nx, params.ny, refine_level);
    fprintf(stderr,"Refinement when relative error is >%Lg%% or absolute error is >%Lg%% of estimated total flux, for pixel brightness > %Lf%% of average\n",
            params.refine_rel*100, params.refine_abs*100, params.refine_cut*100);
  } else {
    fprintf(stderr,"Resolution: %dx%d\n", params.nx, params.ny);
  }

  REAL *taus = calloc(nx*ny, sizeof(*taus));
  REAL *image = calloc(nx*ny, sizeof(*image));
  REAL *imageS = NULL;
  if (taus == NULL || image == NULL) {
    fprintf(stderr, "Could not allocate image memory!");
    exit(-1);
  }
  if(!params.only_unpolarized) {
    imageS = calloc(nx*ny*NIMG, sizeof(*imageS));
    if (imageS == NULL) {
      fprintf(stderr, "Could not allocate image memory!");
      exit(-1);
    }
  }

  // BASE IMAGE at n_min
  // Allocate it, or use the existing allocation for just 1 level
  size_t initialspacingx = (nx - 1) / (nxmin - 1);
  size_t initialspacingy = (ny - 1) / (nymin - 1);

#if DEBUG
  fprintf(stderr, "Image dimensions: %d %d, memory dimensions %ld %ld, minimum %ld %ld\n",
          params.nx, params.ny, nx, ny, nxmin, nymin);

  fprintf(stderr, "Intial spacing: %ld\n", initialspacingx);
#endif

  int *interp_flag = calloc(nx * ny, sizeof(*interp_flag));
  REAL *prelimarray = NULL;
  if (interp_flag == NULL) {
    fprintf(stderr, "Could not allocate adaptive memory!\n");
    exit(-1);
  }
  if (refine_level > 1) {
    prelimarray = calloc(nxmin * nymin, sizeof(*prelimarray));
    if (prelimarray == NULL) {
      fprintf(stderr, "Could not allocate adaptive memory!\n");
      exit(-1);
    }
  }

  // Get the "base image" -- if not adaptively refining, this is just the whole image
  // Note the average of the interpolated image != the average of calculated pixels, though they're close.
#pragma omp parallel for schedule(dynamic,1) collapse(2)
  for (size_t i = 0; i < nx; i += initialspacingx) {
    for (size_t j = 0; j < ny; j += initialspacingy) {

      if (j==0) fprintf(stderr, "%ld ", i);
      size_t thislocation = i / initialspacingy * nymin + j / initialspacingx;
      REAL Intensity = 0;
      REAL Is = 0, Qs = 0, Us = 0, Vs = 0;
      REAL Tau = 0, tauF = 0;

      get_pixel(i, j, params.nx, params.ny, Xcam, params, fovx, fovy, freq,
                  params.only_unpolarized, scale, &Intensity, &Is, &Qs, &Us, &Vs,
                  &Tau, &tauF);

      save_pixel(image, imageS, taus, i, j, nx, ny, params.only_unpolarized, Intensity, Is, Qs, Us,
                  Vs, freqcgs, Tau, tauF);

      interp_flag[i*ny+j] = 0;

      if (refine_level > 1) {
        // computes a total interpolated flux from first pass
        // adds the total number of pixels adjacent to this one
        // assumes nx=ny and nx_min=ny_min
        
        if (i % (nx - 1) != 0 && j % (ny - 1) != 0) {
          // middle case
          prelimarray[thislocation] = Intensity * initialspacingx * initialspacingx;
        } else if ((i == 0 && j % (ny - 1) != 0) || (i % (nx - 1) != 0 && j == 0)) {
          // bottom or left vertical edge
          prelimarray[thislocation] = Intensity * (initialspacingx / 2 + 1) * initialspacingx;
        } else if (i % (nx - 1) != 0 || j % (ny - 1) != 0) {
          // top or right vertical edge
          prelimarray[thislocation] = Intensity * (initialspacingx / 2) * initialspacingx;
        } else {
          // corner case
          if (i == 0 && j == 0) {
            // bottom left corner
            prelimarray[thislocation] = Intensity * (initialspacingx / 2 + 1) * (initialspacingx / 2 + 1);
          } else if (i == 0 || j == 0) {
            // bottom right or top left
            prelimarray[thislocation] = Intensity * (initialspacingx / 2 + 1) * initialspacingx / 2;
          } else {
            // top right
            prelimarray[thislocation] = Intensity * (initialspacingx * initialspacingx) / 4;
          }
        }
      }
    }
  }

  // compute estimated flux total and intensity average
  REAL Iavg = 0;
  if (refine_level > 1) {
    REAL interp_tot = 0;
    for (size_t i = 0; i < nxmin * nymin; i++) {
      interp_tot += prelimarray[i];
    }

    // Average intensity per pixel
    Iavg = interp_tot * pow(freqcgs, 3) / (nx * ny);

    // Print calculated total intensity for debug
#if DEBUG
    fprintf(stderr, "\nInitial flux guess: %Lg", Iavg * (params.nx * params.ny) * scale);
#endif
    fprintf(stderr, "\n\n"); // TODO even for non-refined?
  }

  for (int refined_level = 1; refined_level < refine_level; refined_level++) {
    size_t newspacingx = initialspacingx / pow(2, refined_level);
    size_t newspacingy = initialspacingy / pow(2, refined_level);
    fprintf(stderr, "Refining level %d of %d, spacing %ld,%ld\n", refined_level+1, refine_level, newspacingx, newspacingy);

#pragma omp parallel for schedule(dynamic,1) collapse(2)
    for (size_t i = 0; i < nx; i += newspacingx) {
      for (size_t j = 0; j < ny; j += newspacingy) {

        if (j == 0) fprintf(stderr, "%ld ", i);

        REAL Intensity = 0;
        REAL Is = 0, Qs = 0, Us = 0, Vs = 0;
        REAL Tau = 0, tauF = 0;

        size_t previousspacingx = newspacingx * 2;
        size_t previousspacingy = newspacingy * 2;

        REAL I1, I2, I3, I4, err_abs, err_rel;
        if (i % previousspacingx == 0 && j % previousspacingy == 0) {
          // pixel has already been ray-traced
          continue;
        } else if (i % previousspacingx == 0 && j % previousspacingy != 0) {
          // pixel lies on pre-existing column

          I1 = image[i*ny+j-newspacingy]; // below
          I2 = image[i*ny+j+newspacingy]; // above
          err_abs = (I2 - I1) / 2 / Iavg;
          err_rel = (I2 - I1) / 2 / I1;

          if ((fabs(err_abs) > params.refine_abs && // could be changed to || if wanted
              fabs(err_rel) > params.refine_rel)
              && fabs(I1) > params.refine_cut) {
            // ray trace (tolerances exceeded)

            get_pixel(i, j, params.nx, params.ny, Xcam, params, fovx, fovy, freq,
                        params.only_unpolarized, scale, &Intensity, &Is, &Qs, &Us,
                        &Vs, &Tau, &tauF);

            save_pixel(image, imageS, taus, i, j, nx, ny, params.only_unpolarized, Intensity, Is,
                        Qs, Us, Vs, freqcgs, Tau, tauF);

            interp_flag[i*ny+j] = 0;
                
            } else {
            // interpolate
            interp_flag[i*ny+j] = 1;
            if (params.nearest_neighbor) {
              // nearest interpolation
              save_pixelTransfer(image, imageS, taus, i, j - newspacingy, i,
                                  j, nx, ny, params.only_unpolarized);
              // fills in with the nearest neighbor (choosing one side)
            } else {
              // linear
              lininterp2(image, imageS, taus, i, j - newspacingy, i,
                          j + newspacingy, i, j, nx, ny, params.only_unpolarized);

            }
          }
        } else if (i % previousspacingx != 0 && j % previousspacingy == 0) {
          // pixel lies on pre-existing row
          I1 = image[(i-newspacingx)*ny+j]; //left
          I2 = image[(i+newspacingx)*ny+j]; //right
          err_abs = (I2 - I1) / 2 / Iavg;
          err_rel = (I2 - I1) / 2 / I1;

          if ((fabs(err_abs) > params.refine_abs && //could be changed back to || if wanted
              fabs(err_rel) > params.refine_rel)
              && fabs(I1) > params.refine_cut) {
            // ray trace (tolerances exceeded)

            get_pixel(i, j, params.nx, params.ny, Xcam, params, fovx, fovy, freq,
                        params.only_unpolarized, scale, &Intensity, &Is, &Qs, &Us,
                        &Vs, &Tau, &tauF);

            save_pixel(image, imageS, taus, i, j, nx, ny, params.only_unpolarized, Intensity, Is,
                        Qs, Us, Vs, freqcgs, Tau, tauF);

            interp_flag[i*ny+j] = 0;

          } else {
            // interpolate
            interp_flag[i*ny+j] = 1;
            if (params.nearest_neighbor) {
              // nearest interpolation
              save_pixelTransfer(image, imageS, taus, i - newspacingx, j, i,
                                  j, nx, ny, params.only_unpolarized);
              // fills in with the nearest neighbor (choosing one side)
            } else {
              // linear
              lininterp2(image, imageS, taus, i - newspacingx, j,
                          i + newspacingx, j, i, j, nx, ny, params.only_unpolarized);
            }
          }
        } else {
          // pixel lies equidistant from four corners
          I1 = image[(i-newspacingx)*ny+j-newspacingy]; // bottom left
          I2 = image[(i+newspacingx)*ny+j-newspacingy]; // bottom right
          I3 = image[(i-newspacingx)*ny+j+newspacingy]; // upper left
          I4 = image[(i+newspacingx)*ny+(j+newspacingy)]; // upper right

          // Refinement criterion thanks to Zack Gelles: absolute & relative error of
          // central corner under nearest-neighbor, estimated by Taylor expanding at lower-left pixel
          // Make sure absolute error is in Jy/muas^2

          REAL err_abs = ((I2 + I3) / 2 - I1) / Iavg;
          REAL err_rel = (I2 + I3) / (2 * I1) - 1.;

          if ((fabs (err_abs) > params.refine_abs && //could be changed to && if wanted
              fabs (err_rel) > params.refine_rel)
              && fabs (I1) > params.refine_cut) {

            // ray trace (tolerances exceeded)

            get_pixel(i, j, params.nx, params.ny, Xcam, params, fovx, fovy, freq,
                        params.only_unpolarized, scale, &Intensity, &Is, &Qs, &Us,
                        &Vs, &Tau, &tauF);

            save_pixel(image, imageS, taus, i, j, nx, ny, params.only_unpolarized, Intensity, Is,
                        Qs, Us, Vs, freqcgs, Tau, tauF);

            interp_flag[i*ny+j] = 0;

          } else {
            // interpolate
            interp_flag[i*ny+j] = 1;
            if (params.nearest_neighbor) {
              // nearest interpolation
              save_pixelTransfer(image, imageS, taus, i - newspacingx,
                                  j - newspacingy, i, j, nx, ny, params.only_unpolarized);
              // fills in with the nearest neighbor (choosing one side)
            } else {
              // linear
              lininterp4(image, imageS, taus, i - newspacingx,
                          j - newspacingy, i + newspacingx, j - newspacingy,
                          i - newspacingx, j + newspacingy, i + newspacingx,
                          j + newspacingy, i, j, nx, ny, params.only_unpolarized);
            }
          }
        }
      }
    }

    // Print how many pixels were interpolated
    size_t total_interpolated = 0;
#pragma omp parallel for collapse(2) reduction(+:total_interpolated)
    for (size_t i = 0; i < nx; i++) {
      for (size_t j = 0; j < ny; j++) {
        total_interpolated += interp_flag[i*ny+j];
      }
    }
    // Report interpolation stats vs the number of computed, relevant pixels, not the total
    size_t nx_level = params.nx/newspacingx;
    size_t ny_level = params.ny/newspacingy;
    fprintf(stderr, "\n%ld of %ld (%Lf%%) of computed pixels at %ldx%ld were interpolated\n\n",
            total_interpolated, nx_level * ny_level, ((REAL) total_interpolated) / (nx_level * ny_level) * 100, nx_level, ny_level);
  }

  // TODO print only for "real" pixels
  print_image_stats(image, imageS, nx, ny, params, scale);

  // don't dump if we've been asked to quench output. useful for batch jobs
  // like when fitting light curve fluxes
  if (!params.quench_output) {
    // dump result. if specified, also output ppm image
    dump(image, imageS, taus, params.outf, scale, Xcam, fovx, fovy, nx, ny, &params, params.only_unpolarized);
    if (params.add_ppm) {
      // TODO respect filename from params?
      make_ppm(image, freq, nx, ny, "ipole_lfnu.ppm");
    }
  }

  time = omp_get_wtime() - time;
  printf("Total wallclock time: %Lg s\n\n", time);

  return 0;
}

// TODO Move these?
void get_pixel(size_t i, size_t j, int nx, int ny, REAL Xcam[NDIM], Params params,
               REAL fovx, REAL fovy, REAL freq, int only_intensity, REAL scale,
               REAL *Intensity, REAL *Is, REAL *Qs, REAL *Us, REAL *Vs,
               REAL *Tau, REAL *tauF)
{
  REAL X[NDIM], Kcon[NDIM];
  _Complex REAL N_coord[NDIM][NDIM];

  // Integrate backward to find geodesic trajectory
  init_XK(i,j, params.nx, params.ny, Xcam, params, fovx, fovy, X, Kcon);
  struct of_traj *traj = calloc(params.maxnstep, sizeof(struct of_traj));
#if !INTEGRATOR_TEST
  MULOOP Kcon[mu] *= freq;
#endif
  int nstep = trace_geodesic(X, Kcon, traj, params.eps, params.maxnstep);
  if (nstep >= params.maxnstep-1) {
    // You almost certainly don't want to continue if this happens
    fprintf(stderr, "\nMaxNStep exceeded in pixel %ld %ld!\n", i, j);
    exit(-10);
  }

  // Integrate emission forward along trajectory
  int oddflag = integrate_emission(traj, nstep, Intensity, Tau, tauF, N_coord, &params);

  if (!only_intensity) {
    project_N(X, Kcon, N_coord, Is, Qs, Us, Vs, params.rotcam);
  }

  // Catch anything with bad tetrads, anything we manually specify, and signficantly negative pixels
  if (oddflag != 0 || (i == DIAG_PX_I && j == DIAG_PX_J) || *Is < -1.e-10) {
    fprintf(stderr, "\nOddity in pixel %ld %ld (flag %d):\n", i, j, oddflag);
    //print_vector("Starting X", X);
    print_vector("Starting Kcon", Kcon);
    fprintf(stderr, "nstep: %d\n", nstep);
    fprintf(stderr, "Final Stokes parameters: [%Lg %Lg %Lg %Lg]\n", *Is, *Qs, *Us, *Vs);
  }

  // Record values along the geodesic if requested
  // TODO this is most likely not compatible with adaptive mode
  if (params.trace) {
    size_t stride = params.trace_stride;
    if (params.trace_i < 0 || params.trace_j < 0) { // If no single point is specified
      if (i % stride == 0 && j % stride == 0) { // Save every stride pixels
#pragma omp critical
        dump_var_along(i/stride, j/stride, nstep, traj, params.nx/stride, params.ny/stride, scale, Xcam, fovx, fovy, &params);
      }
    } else {
      if (i == params.trace_i && j == params.trace_j) { // Save just the one
#pragma omp critical
        {
          dump_var_along(0, 0, nstep, traj, 1, 1, scale, Xcam, fovx, fovy, &params);
        }
      }
    }
  }

  free(traj);
}

void save_pixel(REAL *image, REAL *imageS, REAL *taus, size_t i, size_t j, size_t nx, size_t ny, int only_intensity,
                REAL Intensity, REAL Is, REAL Qs, REAL Us, REAL Vs,
                REAL freqcgs, REAL Tau, REAL tauF)
{
  // deposit the intensity and Stokes parameter in pixel
  image[i*ny+j] = Intensity * pow(freqcgs, 3);
  taus[i*ny+j] = Tau;

  if (!only_intensity) {
    imageS[(i*ny+j)*NIMG+0] = Is * pow(freqcgs, 3);
    if (params.qu_conv == 0) {
      imageS[(i*ny+j)*NIMG+1] = -Qs * pow(freqcgs, 3);
      imageS[(i*ny+j)*NIMG+2] = -Us * pow(freqcgs, 3);
    } else {
      imageS[(i*ny+j)*NIMG+1] = Qs * pow(freqcgs, 3);
      imageS[(i*ny+j)*NIMG+2] = Us * pow(freqcgs, 3);
    }
    imageS[(i*ny+j)*NIMG+3] = Vs * pow(freqcgs, 3);
    imageS[(i*ny+j)*NIMG+4] = tauF;

    if (isnan(imageS[(i*ny+j)*NIMG+0])) {
      fprintf(stderr, "NaN in image! Exiting.\n");
      exit(-1);
    }
  }
}

void save_pixelTransfer(REAL *image, REAL *imageS, REAL *taus, size_t iold, size_t jold,
                        size_t inew, size_t jnew, size_t nx, size_t ny, int only_intensity)
{
  // deposit the intensity and Stokes parameter in pixel
    REAL Intensity=image[iold*ny+jold];
    image[inew*ny+jnew]=Intensity;

    if(!only_intensity) {
      REAL Tau=taus[iold*ny+jold];
      REAL Is=imageS[(iold*ny+jold)*NIMG+0];
      REAL Qs=imageS[(iold*ny+jold)*NIMG+1];
      REAL Us=imageS[(iold*ny+jold)*NIMG+2];
      REAL Vs=imageS[(iold*ny+jold)*NIMG+3];
      REAL tauF=imageS[(iold*ny+jold)*NIMG+4];

      taus[inew*ny+jnew] = Tau;
      imageS[(inew*ny+jnew)*NIMG+0] = Is;
      imageS[(inew*ny+jnew)*NIMG+1] = Qs;
      imageS[(inew*ny+jnew)*NIMG+2] = Us;
      imageS[(inew*ny+jnew)*NIMG+3] = Vs;
      imageS[(inew*ny+jnew)*NIMG+4] = tauF;

      if (isnan(imageS[(iold*ny+jold)*NIMG+0])) {
        fprintf(stderr, "NaN in image! Exiting.\n");
        exit(-1);
      }
    }
}

void lininterp2(REAL *image, REAL *imageS, REAL *taus, size_t i1, size_t j1,
                size_t i2, size_t j2, size_t inew, size_t jnew, size_t nx, size_t ny, int only_intensity)
{
  // deposit the intensity and Stokes parameter in pixel
  REAL Intensity = .5 * (image[i1*ny+j1] + image[i2*ny+j2]);
  image[inew*ny+jnew] = Intensity;
  REAL Tau=.5*(taus[i1*ny+j1]+taus[i2*ny+j2]);
  taus[inew*ny+jnew] = Tau;

  if(!only_intensity) {
    REAL Is=.5*(imageS[(i1*ny+j1)*NIMG+0]+imageS[(i2*ny+j2)*NIMG+0]);
    REAL Qs=.5*(imageS[(i1*ny+j1)*NIMG+1]+imageS[(i2*ny+j2)*NIMG+1]);
    REAL Us=.5*(imageS[(i1*ny+j1)*NIMG+2]+imageS[(i2*ny+j2)*NIMG+2]);
    REAL Vs=.5*(imageS[(i1*ny+j1)*NIMG+3]+imageS[(i2*ny+j2)*NIMG+3]);
    REAL tauF=.5*(imageS[(i1*ny+j1)*NIMG+4]+imageS[(i2*ny+j2)*NIMG+4]);

    imageS[(inew*ny+jnew)*NIMG+0] = Is;
    imageS[(inew*ny+jnew)*NIMG+1] = Qs;
    imageS[(inew*ny+jnew)*NIMG+2] = Us;
    imageS[(inew*ny+jnew)*NIMG+3] = Vs;
    imageS[(inew*ny+jnew)*NIMG+4] = tauF;

    if (isnan(imageS[(i1*ny+j1)*NIMG+0])||isnan(imageS[(i2*ny+j2)*NIMG+0])) {
      fprintf(stderr, "NaN in image! Exiting.\n");
      exit(-1);
    }
  }
}


void lininterp4(REAL *image, REAL *imageS, REAL *taus, size_t i1, size_t j1,
                size_t i2,size_t j2, size_t i3, size_t j3, size_t i4, size_t j4, size_t inew, size_t jnew,
                size_t nx, size_t ny, int only_intensity)
{
  // deposit the intensity and Stokes parameter in pixel
  REAL Intensity = .25 * (image[i1*ny+j1] + image[i2*ny+j2] + image[i3*ny+j3] + image[i4*ny+j4]);
  image[inew*ny+jnew] = Intensity;
  REAL Tau=.25*(taus[i1*ny+j1]+taus[i2*ny+j2]+taus[i3*ny+j3]+taus[i4*ny+j4]);
  taus[inew*ny+jnew] = Tau;

  if(!only_intensity) {
    REAL Is=.25*(imageS[(i1*ny+j1)*NIMG+0]+imageS[(i2*ny+j2)*NIMG+0]+imageS[(i3*ny+j3)*NIMG+0]+imageS[(i4*ny+j4)*NIMG+0]);
    REAL Qs=.25*(imageS[(i1*ny+j1)*NIMG+1]+imageS[(i2*ny+j2)*NIMG+1]+imageS[(i3*ny+j3)*NIMG+1]+imageS[(i4*ny+j4)*NIMG+1]);
    REAL Us=.25*(imageS[(i1*ny+j1)*NIMG+2]+imageS[(i2*ny+j2)*NIMG+2]+imageS[(i3*ny+j3)*NIMG+2]+imageS[(i4*ny+j4)*NIMG+2]);
    REAL Vs=.25*(imageS[(i1*ny+j1)*NIMG+3]+imageS[(i2*ny+j2)*NIMG+3]+imageS[(i3*ny+j3)*NIMG+3]+imageS[(i4*ny+j4)*NIMG+3]);
    REAL tauF=.25*(imageS[(i1*ny+j1)*NIMG+4]+imageS[(i2*ny+j2)*NIMG+4]+imageS[(i3*ny+j3)*NIMG+4]+imageS[(i4*ny+j4)*NIMG+4]);

    imageS[(inew*ny+jnew)*NIMG+0] = Is;
    imageS[(inew*ny+jnew)*NIMG+1] = Qs;
    imageS[(inew*ny+jnew)*NIMG+2] = Us;
    imageS[(inew*ny+jnew)*NIMG+3] = Vs;
    imageS[(inew*ny+jnew)*NIMG+4] = tauF;

    if (isnan(imageS[(i1*ny+j1)*NIMG+0])||isnan(imageS[(i2*ny+j2)*NIMG+0])||isnan(imageS[(i3*ny+j3)*NIMG+0])||isnan(imageS[(i4*ny+j4)*NIMG+0])) {
      fprintf(stderr, "NaN in image! Exiting.\n");
      exit(-1);
    }
  }
}

void print_image_stats(REAL *image, REAL *imageS, size_t nx, size_t ny, Params params, REAL scale)
{
  REAL Ftot = 0.;
  REAL Ftot_unpol = 0.;
  REAL Imax = 0.0;
  REAL Iavg = 0.0;
  REAL Qtot = 0.;
  REAL Utot = 0.;
  REAL Vtot = 0.;
  size_t imax = 0;
  size_t jmax = 0;
  for (size_t i = 0; i < params.nx; i++) {
    for (size_t j = 0; j < params.ny; j++) {
      Ftot_unpol += image[i*ny+j]*scale;

      if (!params.only_unpolarized) {
        Ftot += imageS[(i*ny+j)*NIMG+0] * scale;
        Iavg += imageS[(i*ny+j)*NIMG+0];
        Qtot += imageS[(i*ny+j)*NIMG+1] * scale;
        Utot += imageS[(i*ny+j)*NIMG+2] * scale;
        Vtot += imageS[(i*ny+j)*NIMG+3] * scale;
        if (imageS[(i*ny+j)*NIMG+0] > Imax) {
          imax = i;
          jmax = j;
          Imax = imageS[(i*ny+j)*NIMG+0];
        }
      }
    }
  }

  // output normal flux quantities
  fprintf(stderr, "\nscale = %Le\n", scale);
  fprintf(stderr, "imax=%ld jmax=%ld Imax=%Lg Iavg=%Lg\n", imax, jmax, Imax, Iavg/(params.nx*params.ny));
  fprintf(stderr, "freq: %Lg Ftot: %Lg Jy (%Lg Jy unpol xfer) scale=%Lg\n", params.freqcgs, Ftot, Ftot_unpol, scale);
  fprintf(stderr, "nuLnu = %Lg erg/s\n", 4.*M_PI*Ftot * params.dsource * params.dsource * JY * params.freqcgs);

  // output polarized transport information
  REAL LPfrac = 100.*sqrt(Qtot*Qtot+Utot*Utot)/Ftot;
  REAL CPfrac = 100.*Vtot/Ftot;
  fprintf(stderr, "I,Q,U,V [Jy]: %Lg %Lg %Lg %Lg\n", Ftot, Qtot, Utot, Vtot);
  fprintf(stderr, "LP,CP [%%]: %Lg %Lg\n", LPfrac, CPfrac);
}
