#ifndef _PARAMETERS_TD_H
#define _PARAMETERS_TD_H

struct TitovDemoulinParameters
{
  // plasma beta
  double beta;

  // coronal height scale
  double L_G;

  // Torus winding number
  double N_t;

  // Torus major radius
  double R;

  // Magnetic charge separation distance
  double L;

  // Geometrical factor
  double d;

  // The coronal/prominence temperature ratio
  double Tc2Tp;

  // Normalised magnetic charge corresponding to global equilibrium (Eq. 6)
  double q_mag;

  double omega_0;

  double t_drive;

  double t_ramp;
};

#endif