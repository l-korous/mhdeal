#ifndef _PARAMETERS_CS_H
#define _PARAMETERS_CS_H

struct CSParameters
{
  // plasma beta
  double beta;

  // density
  double rho_0;

  // coronal height scale
  double L_G;

  // Torus winding number
  double N_t;

  // Torus major radius
  double R;

  // Torus minor radius
  double r;

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