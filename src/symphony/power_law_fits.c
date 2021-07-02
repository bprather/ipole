#include "power_law.h"

/*power_law_I: fitting formula to the emissivity, in Stokes I, produced by a
 *             power-law distribution of electrons (without any exponential
 *             cutoff).  Uses eq. 29, 33 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to power-law emissivity polarized in Stokes I
 */
REAL power_law_I(struct parameters * params)
{
  REAL nu_c = get_nu_c(*params);

  REAL prefactor = (params->electron_density*pow(params->electron_charge,2.)
                      *nu_c)/params->speed_light;

  REAL term1 = pow(3., params->power_law_p/2.)*(params->power_law_p-1.)
                 *sin(params->observer_angle);

  REAL term2 = 2.*(params->power_law_p+1.)
                 *(pow(params->gamma_min, 1.-params->power_law_p)
                   -pow(params->gamma_max, 1.-params->power_law_p));

  REAL term3 = tgamma((3.*params->power_law_p-1.)/12.)
                *tgamma((3.*params->power_law_p+19.)/12.);

  REAL term4 = pow(params->nu/(nu_c*sin(params->observer_angle)), 
                     -(params->power_law_p-1.)/2.);

  REAL ans = prefactor*term1/term2*term3*term4;

  return ans;
}

/*power_law_Q: fitting formula to the emissivity, in Stokes Q, produced by a
 *             power-law distribution of electrons (without any exponential
 *             cutoff).  Uses eq. 29, 33 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to power-law emissivity polarized in Stokes Q
 */
REAL power_law_Q(struct parameters * params)
{
  REAL p_term = -(params->power_law_p + 1.)/(params->power_law_p + 7./3.);

  REAL ans = p_term * power_law_I(params);

  return ans;
}

/*power_law_V: fitting formula to the emissivity, in Stokes V, produced by a
 *             power-law distribution of electrons (without any exponential
 *             cutoff).  Uses eq. 29, 33 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to power-law emissivity polarized in Stokes V
 */
REAL power_law_V(struct parameters * params)
{
  REAL nu_c = get_nu_c(*params);

  REAL term1 = -(171./250.)*pow(params->power_law_p, 49./100.);

  REAL term2 = 1./tan(params->observer_angle) 
                 * pow(params->nu
                       /(3.*nu_c*sin(params->observer_angle)), -1./2.);
  
  REAL ans = term1*term2*power_law_I(params);
 
  /*NOTE: Sign corrected; the sign in Leung et al. (2011)
    and Pandya et al. (2016) for Stokes V transfer coefficients
    does not follow the convention the papers describe (IEEE/IAU);
    the sign has been corrected here.*/ 
  return -ans;
}

/*power_law_I_abs: fitting formula to the absorptivity, in Stokes I, from
 *                 by a power-law distribution of electrons (without any 
 *                 exponential cutoff).  Uses eq. 30, 34 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to power-law absorptivity polarized in Stokes I
 */
REAL power_law_I_abs(struct parameters * params)
{
  REAL nu_c = get_nu_c(*params);

  REAL prefactor = (params->electron_density*pow(params->electron_charge,2.))
                    /(params->nu*params->mass_electron*params->speed_light);

  REAL term1 = pow(3., (params->power_law_p+1.)/2.)*(params->power_law_p-1.);

  REAL term2 = 4.*(pow(params->gamma_min, 1.-params->power_law_p)
                -pow(params->gamma_max, 1.-params->power_law_p));

  REAL term3 = tgamma((3.*params->power_law_p+2.)/12.)
                *tgamma((3.*params->power_law_p+22.)/12.);

  REAL term4 = pow(params->nu/(nu_c*sin(params->observer_angle)), 
                     -(params->power_law_p+2.)/2.);

  REAL ans = prefactor*term1/term2*term3*term4;

  return ans;
}

/*power_law_Q_abs: fitting formula to the absorptivity, in Stokes Q, from
 *                 by a power-law distribution of electrons (without any 
 *                 exponential cutoff).  Uses eq. 30, 34 of [1].
 *
 *@params: struct of parameters params
 *@returns: fitting formula to power-law absorptivity polarized in Stokes Q
 */
REAL power_law_Q_abs(struct parameters * params)
{
  REAL nu_c = get_nu_c(*params);

  REAL prefactor = (params->electron_density*pow(params->electron_charge,2.))
                    /(params->nu*params->mass_electron*params->speed_light);

  REAL term1 = pow(3., (params->power_law_p+1.)/2.)*(params->power_law_p-1.);

  REAL term2 = 4.*(pow(params->gamma_min, 1.-params->power_law_p)
                -pow(params->gamma_max, 1.-params->power_law_p));

  REAL term3 = tgamma((3.*params->power_law_p+2.)/12.)
                *tgamma((3.*params->power_law_p+22.)/12.);

  REAL term4 = pow(params->nu/(nu_c*sin(params->observer_angle)), 
                     -(params->power_law_p+2.)/2.);

  REAL term5 = -pow((17./500.)*params->power_law_p - 43./1250., 43./500.);

  REAL ans = prefactor*term1/term2*term3*term4*term5;

  return ans;
}

/*power_law_V_abs: fitting formula to the absorptivity, in Stokes V, from
 *                 by a power-law distribution of electrons (without any 
 *                 exponential cutoff).  Uses eq. 30, 34 of [1], but
 *                 is modified slightly for increased accuracy at the
 *                 cost of a more complicated function.
 *
 *@params: struct of parameters params
 *@returns: fitting formula to power-law absorptivity polarized in Stokes V
 */
REAL power_law_V_abs(struct parameters * params)
{
  REAL nu_c = get_nu_c(*params);

  REAL prefactor = (params->electron_density*pow(params->electron_charge,2.))
                    /(params->nu*params->mass_electron*params->speed_light);

  REAL term1 = pow(3., (params->power_law_p+1.)/2.)*(params->power_law_p-1.);

  REAL term2 = 4.*(pow(params->gamma_min, 1.-params->power_law_p)
                -pow(params->gamma_max, 1.-params->power_law_p));

  REAL term3 = tgamma((3.*params->power_law_p+2.)/12.)
                *tgamma((3.*params->power_law_p+22.)/12.);

  REAL term4 = pow(params->nu/(nu_c*sin(params->observer_angle)), 
                     -(params->power_law_p+2.)/2.);

  REAL term5 = -pow((71./100.)*params->power_law_p+22./625.,197./500.);

  REAL term6 = pow((31./10.)*pow(sin(params->observer_angle),-48./25)-31./10.,
                     64./125.);

  REAL term7 = pow(params->nu/(nu_c*sin(params->observer_angle)), -1./2.);
 
  REAL ans = prefactor*term1/term2*term3*term4*term5*term6*term7;

  /*The Stokes V absorption coefficient changes sign at observer_angle
    equals 90deg, but this formula does not.  This discrepancy is a 
    bug in this formula, and is patched by the term below.*/
  REAL sign_bug_patch = cos(params->observer_angle) / 
                          fabs(cos(params->observer_angle));

  /*NOTE: Sign corrected; the sign in Leung et al. (2011)
    and Pandya et al. (2016) for Stokes V transfer coefficients
    does not follow the convention the papers describe (IEEE/IAU);
    the sign has been corrected here.*/
  return -ans * sign_bug_patch;
}

