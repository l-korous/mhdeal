#include "completeEllipticIntegrals.h"
#include "boundaryConditionTD.h"
#include "equationsMhd.h"

template <int dim>
BoundaryConditionTDWithVortices<dim>::BoundaryConditionTDWithVortices(Parameters<dim>& parameters, TitovDemoulinParameters& td_parameters) :
  BoundaryCondition<EquationsTypeMhd, dim>(parameters), td_parameters(td_parameters)
{
  this->eps = sqrt(1. - (td_parameters.d / td_parameters.R) * (td_parameters.d / td_parameters.R));
  y_1 = this->eps * td_parameters.R;
  y_2 = -this->eps * td_parameters.R;
}

template <int dim>
double BoundaryConditionTDWithVortices<dim>::r_1_bar(double x, double y) const
{
  return sqrt(this->eps * this->eps * ((y - this->y_1) * (y - this->y_1)) + (x * x));
}

template <int dim>
double BoundaryConditionTDWithVortices<dim>::r_2_bar(double x, double y) const
{
  return sqrt(this->eps * this->eps * ((y - this->y_2) * (y - this->y_2)) + (x * x));
}

template <int dim>
double BoundaryConditionTDWithVortices<dim>::omega_1(double x, double y) const
{
  return 0.5 * (1. + tanh(1. - r_1_bar(x, y)));
}

template <int dim>
double BoundaryConditionTDWithVortices<dim>::omega_2(double x, double y) const
{
  return 0.5 * (1. + tanh(1. - r_2_bar(x, y)));
}

template <int dim>
double BoundaryConditionTDWithVortices<dim>::omega(double time) const
{
  return (this->td_parameters.omega_0 / 2.) * (1. + tanh((time - this->td_parameters.t_drive) / this->td_parameters.t_ramp));
}

template <int dim>
void BoundaryConditionTDWithVortices<dim>::bc_vector_value(int boundary_no, const Point<dim> &point, InputVector &result, const InputVector &W_plus, double time) const
{
  // For other than z=0 boundaries, we use do-nothing
  if (point[2] > SMALL)
  {
    for (unsigned int di = 0; di < Equations<EquationsTypeMhd, dim>::n_components; ++di)
      result[di] = W_plus[di];
    return;
  }

  double x = point[0], y = point[1];
  result[0] = W_plus[0];
  result[1] = result[0] * (-omega(time) * this->eps * ((y - y_1) * omega_1(x, y) + (y - y_2) * omega_2(x, y)));
  result[2] = result[0] * (omega(time) * (x / this->eps) * (omega_1(x, y) + omega_2(x, y)));
  result[3] = 0.0;

  // energy density
  result[4] = W_plus[4];

  result[5] = W_plus[5];
  result[6] = W_plus[6];
  result[7] = W_plus[7];
}

template <int dim>
BoundaryConditionTDFree<dim>::BoundaryConditionTDFree(Parameters<dim>& parameters, TitovDemoulinParameters& td_parameters) :
  BoundaryCondition<EquationsTypeMhd, dim>(parameters)
{
}

template <int dim>
void BoundaryConditionTDFree<dim>::bc_vector_value(int boundary_no, const Point<dim> &point, InputVector &result, const InputVector &W_plus, double time) const
{
  for (unsigned int di = 0; di < Equations<EquationsTypeMhd, dim>::n_components; ++di)
    result[di] = W_plus[di];
  return;
}

template <int dim>
BoundaryConditionTDInitialState<dim>::BoundaryConditionTDInitialState(Parameters<dim>& parameters, TitovDemoulinParameters& td_parameters) :
  BoundaryCondition<EquationsTypeMhd, dim>(parameters), td_parameters(td_parameters)
{
  invL_G = 0.0;
  if (td_parameters.L_G > 0.0)
    invL_G = 1.0 / td_parameters.L_G;

  // Submerging of torus main axis in units of R
  d2R = td_parameters.d / td_parameters.R;

  // Distance of magnetic charges from x=0 plane in units of R
  L2R = td_parameters.L / td_parameters.R;

  R2L = td_parameters.R / td_parameters.L;

  //======================= Calculate dependent TD model parameters
  // Normalised magnetic charge corresponding to global equilibrium (Eq. 6)
  q_mag = 0.25 * fabs(td_parameters.N_t) * (log(8.0 * td_parameters.R) - 1.25) * (1 + L2R * L2R) * sqrt(1 + L2R * L2R) / L2R;
  // Should be ? q_mag = 0.25 * fabs(td_parameters.N_t) * (log(8.0 * td_parameters.R) - 1.25) * (1 + R2L * R2L) * sqrt(1 + R2L * R2L) / (L2R * L2R);
  // This would reflect q_{maq} to be |q| from page 10/456 of https://github.com/l-korous/doctoral-thesis/blob/master/_reference/Modeling%20of%20H%CE%B1%20Eruptive%20Events%20Observed%20at%20the%20Solar.pdf

  // Sign of winding: corresponds to sign of I_O in TD paper
  iSgn = (td_parameters.N_t >= 0) ? 1.0 : -1.0;
}

template <int dim>
void BoundaryConditionTDInitialState<dim>::bc_vector_value(int boundary_no, const Point<dim> &p, InputVector &result, const InputVector &W_plus, double time) const
{
  // For other than z=0 boundaries, we use do-nothing
  if (p[2] > SMALL)
  {
    for (unsigned int di = 0; di < Equations<EquationsTypeMhd, dim>::n_components; ++di)
      result[di] = W_plus[di];
    return;
  }

  //========== Calculate the vector potential for I_t-generated toroidal field 
  double xx, yy, zz;
  Point<dim> &ca = this->parameters.corner_a;
  Point<dim> &cb = this->parameters.corner_b;

  xx = p[0] - (ca[0] + cb[0]) * 0.5;
  yy = p[1] - (ca[1] + cb[1]) * 0.5;
  zz = p[2] - ca[2];

  Vector<double> theta0(3);
  double r_maj, r_min;
  const double dd = 1e-8;  // precision of numerical derivatives
  const double idd = 0.5 / dd;
  double df[6][3];
  double P[6][3];
  double dP[3][3];

  // calculate vector potentil in 6 close points (+dx,-dx,+dy,-dy,+dz,-dz)
  for (unsigned int i = 0; i < 6; ++i) {
    df[i][0] = df[i][1] = df[i][2] = 0.0;
    df[i][int(i / 2)] = double(1 - 2 * int(i % 2)) * dd;
  }

  for (unsigned int i = 0; i < 6; ++i) {
    double x = xx + df[i][0];
    double y = yy + df[i][1];
    double z = zz + df[i][2];

    // Distances from torus major and minor axes
    r_maj = sqrt(y * y + (z + d2R * this->td_parameters.R) * (z + d2R * this->td_parameters.R));
    r_min = sqrt(x * x + (r_maj - this->td_parameters.R) * (r_maj - this->td_parameters.R));

    // Unit vector of toroidal coordinate theta
    theta0[0] = 0.0;
    theta0[1] = -(z + d2R * this->td_parameters.R) / r_maj;
    theta0[2] = y / r_maj;

    // Common radial factor for A_tor
    double rFactor = fabs(this->td_parameters.N_t) * sqrt(1.0 / (this->td_parameters.R * r_maj));

    // Argument of elliptical integral
    double kr = 2.0 * sqrt(r_maj * this->td_parameters.R / ((r_maj + this->td_parameters.R) * (r_maj + this->td_parameters.R) + x * x));

    //---- Sew-up internal and external solutions
    if (r_min > 1.0) { //---- external region 
                     // Elliptical integrals
      double Ek, Kk;
      Complete_Elliptic_Integrals_Modulus(kr, Kk, Ek);

      double Ak = ((2.0 - kr * kr) * Kk - 2.0 * Ek) / kr;

      P[i][0] = (rFactor * Ak) * theta0[0];
      P[i][1] = (rFactor * Ak) * theta0[1];
      P[i][2] = (rFactor * Ak) * theta0[2];

    }
    else { //---- inside the torus

           // ka=kr at r_min=1 (=torus surface) 
      double ka = 2.0 * sqrt(r_maj * td_parameters.R / (4.0 * r_maj * this->td_parameters.R + 1.0));

      double Ek, Kk;
      Complete_Elliptic_Integrals_Modulus(ka, Kk, Ek);

      double Ak = ((2.0 - ka * ka) * Kk - 2.0 * Ek) / ka;
      double Ak_prime = ((2.0 - ka * ka) * Ek - 2.0 * (1.0 - ka * ka) * Kk) /
        (ka * ka * (1.0 - ka * ka));

      double cf = (rFactor * (Ak + Ak_prime*(kr - ka)));
      P[i][0] = cf*theta0[0];
      P[i][1] = cf*theta0[1];
      P[i][2] = cf*theta0[2];
    }
  }

  // calculate derivatives of vector potential
  for (unsigned int i = 0; i < 3; ++i) {
    dP[i][0] = (P[2 * i][0] - P[2 * i + 1][0]) * idd;
    dP[i][1] = (P[2 * i][1] - P[2 * i + 1][1]) * idd;
    dP[i][2] = (P[2 * i][2] - P[2 * i + 1][2]) * idd;
  }

  //====================== Calculate the full state field

  double pressure;

  // Distances from torus major and minor axes
  r_maj = sqrt(yy * yy + (zz + d2R * this->td_parameters.R) * (zz + d2R * this->td_parameters.R));
  r_min = sqrt(xx * xx + (r_maj - this->td_parameters.R) * (r_maj - this->td_parameters.R));

  // Unit vector of toroidal coordinate theta
  theta0[0] = 0.0;
  theta0[1] = -(zz + d2R * this->td_parameters.R) / r_maj;
  theta0[2] = yy / r_maj;

  // Radius vectors originating in magnetic charges
  Vector<double> r_plus(3);//(xx-L2R * R,yy,zz+d2R * R);
  r_plus[0] = xx - L2R * this->td_parameters.R;
  r_plus[1] = yy;
  r_plus[2] = zz + d2R * this->td_parameters.R;
  Vector<double> r_minus(3);//(xx+L2R * R,yy,zz+d2R * R);
  r_minus[0] = xx + L2R * this->td_parameters.R;
  r_minus[1] = yy;
  r_minus[2] = zz + d2R * this->td_parameters.R;

  double rp = r_plus.l2_norm();
  double rm = r_minus.l2_norm();

  //========== Calculate the magnetic field in TD equilibrium by parts

  //--- Q-generated part of field

  double cf1 = q_mag / (rp * rp * rp);
  double cf2 = q_mag / (rm * rm * rm);
  Vector<double> B_loc(3);//=cf1*r_plus-cf2*r_minus;
  B_loc.sadd(0.0, cf1, r_plus, -cf2, r_minus);

  // add vector potential part B = curl A
  B_loc[0] += dP[1][2] - dP[2][1];
  B_loc[1] += dP[2][0] - dP[0][2];
  B_loc[2] += dP[0][1] - dP[1][0];

  /*
  barta@asu.cas.cz
  10/02/2012

  With the following density-scaling the velocities are normalised
  to V_A outside the flux-rope; L_G is the coronal height-scale.

  ----------

  The density (and consequent temperature) jump rho~H(r_min-1)
  replaced by smoother rho~tgh(r_min-1) profile (TPCR-like).
  */

  double rho_0 = 0.5 * (1.0 - td_parameters.Tc2Tp) * tanh(densGrad*(r_min - 1.0)) + 0.5 * (1 + td_parameters.Tc2Tp);

  if (r_min > 1.0) { // external region

    result[0] = rho_0 * exp(-zz * invL_G);         // mass density outside

    B_loc.sadd(1.0, iSgn * td_parameters.R / r_maj, theta0);

    pressure = td_parameters.beta * exp(-zz * invL_G);
  }
  else { // inside the torus

    result[0] = rho_0 * exp(-zz * td_parameters.Tc2Tp * invL_G);   // mass density in the loop

    B_loc.sadd(1.0, (iSgn * (sqrt(1.0 + td_parameters.H * (1.0 - r_min * r_min)) + td_parameters.R / r_maj - 1.0)), theta0);

    pressure = td_parameters.beta * exp(-zz * invL_G);
  }

  result[1] = 0.0;
  result[2] = 0.0;                          // momentum density
  result[3] = 0.0;

  // energy density
  result[4] = pressure / (this->parameters.gas_gamma - 1.0) + 0.5 * (B_loc[0] * B_loc[0] + B_loc[1] * B_loc[1] + B_loc[2] * B_loc[2]);

  result[5] = B_loc[0];
  result[6] = B_loc[1];                  // magnetic field
  result[7] = B_loc[2];
}

template class BoundaryConditionTDWithVortices<3>;
template class BoundaryConditionTDFree<3>;
template class BoundaryConditionTDInitialState<3>;
