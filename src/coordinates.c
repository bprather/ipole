// Coordinate-dependent (but not radiation model-dependent) functions for ipole

#include "coordinates.h"

#include "decs.h"
#include "geometry.h"

int use_eKS_internal;
int metric;
REAL a, hslope; // mks
REAL poly_norm, poly_xt, poly_alpha, mks_smooth; // fmks
REAL mks3R0, mks3H0, mks3MY1, mks3MY2, mks3MP0; // mks3
REAL startx[NDIM], stopx[NDIM], dx[NDIM];
REAL R0, Rin, Rout, Rh;

/*
 * Despite the name, this returns r, th coordinates for a KS or BL
 * coordinate system (since they're equal), from a set of "modified"
 * native coordinates X
 * If METRIC_MINKOWSKI is set, returns spherical coordinates instead
 */
void bl_coord(REAL X[NDIM], REAL *r, REAL *th)
{

  if (metric == METRIC_MINKOWSKI) {
      *r = X[1];
      *th = X[2];
      return;
  }

  *r = exp(X[1]);

  if (use_eKS_internal) {
    *th = M_PI * X[2];
  } else {
    REAL y, thG, thJ;
    switch (metric) {
      case METRIC_MKS:
        *th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
        break;
      case METRIC_BHACMKS:
        *th = X[2] + (hslope / 2.) * sin(2. * X[2]);
        break;
      case METRIC_FMKS:
        thG = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
        y = 2 * X[2] - 1.;
        thJ = poly_norm * y
            * (1. + pow(y / poly_xt, poly_alpha) / (poly_alpha + 1.)) + 0.5 * M_PI;
        *th = thG + exp(mks_smooth * (startx[1] - X[1])) * (thJ - thG);
        break;
      case METRIC_MKS3:
        *r = exp(X[1]) + mks3R0;
        *th = (M_PI
            * (1.
                + 1. / tan((mks3H0 * M_PI) / 2.)
                    * tan(
                        mks3H0 * M_PI
                            * (-0.5
                                + (mks3MY1
                                    + (pow(2., mks3MP0) * (-mks3MY1 + mks3MY2))
                                        / pow(exp(X[1]) + mks3R0, mks3MP0))
                                    * (1. - 2. * X[2]) + X[2])))) / 2.;
        break;
    }
  }
}

void bl_to_ks(REAL X[NDIM], REAL ucon_bl[NDIM], REAL ucon_ks[NDIM])
{
  REAL r, th;
  bl_coord(X, &r, &th);

  REAL trans[NDIM][NDIM];
  MUNULOOP
    trans[mu][nu] = delta(mu, nu);

  trans[0][1] = 2. * r / (r * r - 2. * r + a * a);
  trans[3][1] = a / (r * r - 2. * r + a * a);

  MULOOP
    ucon_ks[mu] = 0.;
  MUNULOOP
    ucon_ks[mu] += trans[mu][nu] * ucon_bl[nu];
}

void ks_to_bl(REAL X[NDIM], REAL ucon_ks[NDIM], REAL ucon_bl[NDIM])
{
  REAL r, th;
  bl_coord(X, &r, &th);

  REAL trans[NDIM][NDIM], rev_trans[NDIM][NDIM];
  MUNULOOP
    trans[mu][nu] = delta(mu, nu);

  trans[0][1] = 2. * r / (r * r - 2. * r + a * a);
  trans[3][1] = a / (r * r - 2. * r + a * a);

  invert_matrix(trans, rev_trans);

  MULOOP
    ucon_bl[mu] = 0.;
  MUNULOOP
    ucon_bl[mu] += rev_trans[mu][nu] * ucon_ks[nu];
}

/*
 * returns g_{munu} at location specified by X
 */
void gcov_func(REAL X[NDIM], REAL gcov[NDIM][NDIM])
{
  REAL r, th;
  bl_coord(X, &r, &th);

  if (metric == METRIC_MINKOWSKI) {
      MUNULOOP gcov[mu][nu] = 0;
      gcov[0][0] = -1;
      gcov[1][1] = 1;
      gcov[2][2] = r*r;
      gcov[3][3] = r*r*sin(th)*sin(th);
      return;
  }

  // compute ks metric
  REAL Gcov_ks[NDIM][NDIM];
  gcov_ks(r, th, Gcov_ks);

  // convert from ks metric to mks/mmks
  REAL dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);

  MUNULOOP
  {
    gcov[mu][nu] = 0;
    for (int lam = 0; lam < NDIM; ++lam) {
      for (int kap = 0; kap < NDIM; ++kap) {
        gcov[mu][nu] += Gcov_ks[lam][kap] * dxdX[lam][mu] * dxdX[kap][nu];
      }
    }
  }
}

// compute KS metric at point (r,th) in KS coordinates (cyclic in t, ph)
inline void gcov_ks(REAL r, REAL th, REAL gcov[NDIM][NDIM])
{
  REAL cth = cos(th);
  REAL sth = sin(th);

  REAL s2 = sth * sth;
  REAL rho2 = r * r + a * a * cth * cth;

  MUNULOOP gcov[mu][nu] = 0.;
  // Compute KS metric from KS coordinates (cyclic in t,phi)
  gcov[0][0] = -1. + 2. * r / rho2;
  gcov[0][1] = 2. * r / rho2;
  gcov[0][3] = -2. * a * r * s2 / rho2;

  gcov[1][0] = gcov[0][1];
  gcov[1][1] = 1. + 2. * r / rho2;
  gcov[1][3] = -a * s2 * (1. + 2. * r / rho2);

  gcov[2][2] = rho2;

  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];
  gcov[3][3] = s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2));
}

inline void gcov_bl(REAL r, REAL th, REAL gcov[NDIM][NDIM])
{
  REAL sth, cth, s2, a2, r2, DD, mu;
  sth = fabs(sin(th));
  s2 = sth * sth;
  cth = cos(th);
  a2 = a * a;
  r2 = r * r;
  DD = 1. - 2. / r + a2 / r2;
  mu = 1. + a2 * cth * cth / r2;

  MUNULOOP gcov[mu][nu] = 0.;
  // Compute BL metric from BL coordinates
  gcov[0][0] = -(1. - 2. / (r * mu));
  gcov[0][3] = -2. * a * s2 / (r * mu);
  gcov[3][0] = gcov[0][3];
  gcov[1][1] = mu / DD;
  gcov[2][2] = r2 * mu;
  gcov[3][3] = r2 * sth * sth * (1. + a2 / r2 + 2. * a2 * s2 / (r2 * r * mu));

}

void set_dxdX(REAL X[NDIM], REAL dxdX[NDIM][NDIM])
{
  // Jacobian with respect to KS basis where X is given in
  // non-KS basis
  MUNULOOP
    dxdX[mu][nu] = delta(mu, nu);

  dxdX[1][1] = exp(X[1]); // Overridden by one case below

  if (use_eKS_internal) {
    dxdX[2][2] = M_PI;
  } else {
    switch (metric) {
      case METRIC_MKS:
        dxdX[2][2] = M_PI + (1 - hslope) * M_PI * cos(2. * M_PI * X[2]);
        break;
      case METRIC_BHACMKS:
        dxdX[2][2] = 1 + hslope * cos(2. * X[2]);
        break;
      case METRIC_FMKS:
        dxdX[2][1] = -exp(mks_smooth * (startx[1] - X[1])) * mks_smooth
            * (
            M_PI / 2. -
            M_PI * X[2]
                + poly_norm * (2. * X[2] - 1.)
                    * (1
                        + (pow((-1. + 2 * X[2]) / poly_xt, poly_alpha))
                            / (1 + poly_alpha))
                - 1. / 2. * (1. - hslope) * sin(2. * M_PI * X[2]));
        dxdX[2][2] = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2])
            + exp(mks_smooth * (startx[1] - X[1]))
                * (-M_PI
                    + 2. * poly_norm
                        * (1.
                            + pow((2. * X[2] - 1.) / poly_xt, poly_alpha)
                                / (poly_alpha + 1.))
                    + (2. * poly_alpha * poly_norm * (2. * X[2] - 1.)
                        * pow((2. * X[2] - 1.) / poly_xt, poly_alpha - 1.))
                        / ((1. + poly_alpha) * poly_xt)
                    - (1. - hslope) * M_PI * cos(2. * M_PI * X[2]));
        break;
      case METRIC_MKS3:
        dxdX[2][1] = -(pow(2., -1. + mks3MP0) * exp(X[1]) * mks3H0 * mks3MP0
            * (mks3MY1 - mks3MY2) * pow(M_PI, 2)
            * pow(exp(X[1]) + mks3R0, -1 - mks3MP0) * (-1 + 2 * X[2]) * 1.
            / tan((mks3H0 * M_PI) / 2.)
            * pow(
                1.
                    / cos(
                        mks3H0 * M_PI
                            * (-0.5
                                + (mks3MY1
                                    + (pow(2, mks3MP0) * (-mks3MY1 + mks3MY2))
                                        / pow(exp(X[1]) + mks3R0, mks3MP0))
                                    * (1 - 2 * X[2]) + X[2])),
                2));
        dxdX[2][2] = (mks3H0 * pow(M_PI, 2)
            * (1
                - 2
                    * (mks3MY1
                        + (pow(2, mks3MP0) * (-mks3MY1 + mks3MY2))
                            / pow(exp(X[1]) + mks3R0, mks3MP0))) * 1.
            / tan((mks3H0 * M_PI) / 2.)
            * pow(
                1.
                    / cos(
                        mks3H0 * M_PI
                            * (-0.5
                                + (mks3MY1
                                    + (pow(2, mks3MP0) * (-mks3MY1 + mks3MY2))
                                        / pow(exp(X[1]) + mks3R0, mks3MP0))
                                    * (1 - 2 * X[2]) + X[2])),
                2)) / 2.;
        break;
      case METRIC_MINKOWSKI:
        // Blank transform: just override L_11
        dxdX[1][1] = 1.;
    }
  }
}

void set_dXdx(REAL X[NDIM], REAL dXdx[NDIM][NDIM]) {
  REAL dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);
  invert_matrix(dxdX, dXdx);
}

void vec_to_ks(REAL X[NDIM], REAL v_nat[NDIM], REAL v_ks[NDIM]) {
  REAL dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);

  MULOOP v_ks[mu] = 0.;
  MUNULOOP v_ks[mu] += dxdX[mu][nu] * v_nat[nu];
}

void vec_from_ks(REAL X[NDIM], REAL v_ks[NDIM], REAL v_nat[NDIM]) {
  REAL dXdx[NDIM][NDIM];
  set_dXdx(X, dXdx);

  MULOOP v_nat[mu] = 0.;
  MUNULOOP v_nat[mu] += dXdx[mu][nu] * v_ks[nu];
}

/*
 * Translate the input camera angles into a canonical Xcam in native coordinates
 */
void native_coord(REAL r, REAL th, REAL phi, REAL X[NDIM]) {
  if (metric == METRIC_MINKOWSKI) {
    X[0] = 1; X[1] = r; X[2] = th/180*M_PI; X[3] = phi/180*M_PI;
  } else {
    REAL x[NDIM] = {0., r, th/180.*M_PI, phi/180.*M_PI};
    X[0] = 0.0;
    X[1] = log(r);
    X[2] = root_find(x);
    X[3] = phi/180.*M_PI;
  }
}

/**
 * Root-find the camera theta in native coordinates
 * 
 * TODO switch this to native_coord func above...
 */
REAL root_find(REAL X[NDIM])
{
  REAL th = X[2];
  REAL tha, thb, thc;

  REAL Xa[NDIM], Xb[NDIM], Xc[NDIM];
  Xa[1] = log(X[1]);
  Xa[3] = X[3];
  Xb[1] = Xa[1];
  Xb[3] = Xa[3];
  Xc[1] = Xa[1];
  Xc[3] = Xa[3];

  if (X[2] < M_PI / 2.) {
    Xa[2] = startx[2];
    Xb[2] = (stopx[2] - startx[2])/2 + SMALL;
  } else {
    Xa[2] = (stopx[2] - startx[2])/2 - SMALL;
    Xb[2] = stopx[2];
  }

  REAL tol = 1.e-9;
  tha = theta_func(Xa);
  thb = theta_func(Xb);

  // check limits first
  if (fabs(tha-th) < tol) {
    return Xa[2];
  } else if (fabs(thb-th) < tol) {
    return Xb[2];
  }

  // bisect for a bit
  for (int i = 0; i < 1000; i++) {
    Xc[2] = 0.5 * (Xa[2] + Xb[2]);
    thc = theta_func(Xc);

    if ((thc - th) * (thb - th) < 0.)
      Xa[2] = Xc[2];
    else
      Xb[2] = Xc[2];

    REAL err = theta_func(Xc) - th;
    if (fabs(err) < tol)
      break;
  }

  return Xc[2];
}
