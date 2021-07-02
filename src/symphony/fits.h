#ifndef SYMPHONY_FITS_H_
#define SYMPHONY_FITS_H_

#include "params.h"

/* Maxwell-Juettner emissivity fits */
REAL maxwell_juettner_I(struct parameters * params);
REAL maxwell_juettner_Q(struct parameters * params);
REAL maxwell_juettner_V(struct parameters * params);

/* Maxwell-Juettner absorptivity fits */
REAL maxwell_juettner_I_abs(struct parameters * params);
REAL maxwell_juettner_Q_abs(struct parameters * params);
REAL maxwell_juettner_V_abs(struct parameters * params);

/* Maxwell-Juettner Faraday rotation fits */
REAL maxwell_juettner_rho_Q(struct parameters * params);
REAL maxwell_juettner_rho_V(struct parameters * params);

/* Power-law emissivity fits */
REAL power_law_I(struct parameters * params);
REAL power_law_Q(struct parameters * params);
REAL power_law_V(struct parameters * params);

/* Power-law absorptivity fits */
REAL power_law_I_abs(struct parameters * params);
REAL power_law_Q_abs(struct parameters * params);
REAL power_law_V_abs(struct parameters * params);

/* Kappa emissivity fits */
REAL kappa_I(struct parameters * params);
REAL kappa_Q(struct parameters * params);
REAL kappa_V(struct parameters * params);

/* Kappa absorptivity fits */
REAL kappa_I_abs(struct parameters * params);
REAL kappa_Q_abs(struct parameters * params);
REAL kappa_V_abs(struct parameters * params);

REAL check_for_errors(struct parameters * params);
REAL j_nu_fit(REAL nu,
                REAL magnetic_field,
                REAL electron_density,
                REAL observer_angle,
                int distribution,
                int polarization,
                REAL theta_e,
                REAL power_law_p,
                REAL gamma_min,
                REAL gamma_max,
                REAL gamma_cutoff,
                REAL kappa,
                REAL kappa_width
                );
REAL alpha_nu_fit(REAL nu,
                    REAL magnetic_field,
                    REAL electron_density,
                    REAL observer_angle,
                    int distribution,
                    int polarization,
                    REAL theta_e,
                    REAL power_law_p,
                    REAL gamma_min,
                    REAL gamma_max,
                    REAL gamma_cutoff,
                    REAL kappa,
                    REAL kappa_width
                    );

REAL rho_nu_fit(REAL nu,
                  REAL magnetic_field,
                  REAL electron_density,
                  REAL observer_angle,
                  int distribution,
                  int polarization,
                  REAL theta_e,
                  REAL power_law_p,
                  REAL gamma_min,
                  REAL gamma_max,
                  REAL gamma_cutoff,
                  REAL kappa,
                  REAL kappa_width
		  );
#endif /* SYMPHONY_FITS_H_ */

