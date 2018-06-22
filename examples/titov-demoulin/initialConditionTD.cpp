#include "completeEllipticIntegrals.h"
#include "initialConditionTD.h"
#include "equationsMhd.h"

template <EquationsType equationsType, int dim>
InitialConditionTitovDemoulin<equationsType, dim>::InitialConditionTitovDemoulin(Parameters<dim>& parameters, TitovDemoulinParameters& td_parameters) :
  InitialCondition<equationsType, dim>(parameters), td_parameters(td_parameters)
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
  // q_mag = 0.25 * fabs(td_parameters.N_t) * (log(8.0 * td_parameters.R) - 1.25) * (1 + L2R * L2R) * sqrt(1 + L2R * L2R) / L2R;
  q_mag = 0.25 * fabs(td_parameters.N_t) * (log(8.0 * td_parameters.R) - 1.25) * (1 + R2L * R2L) * sqrt(1 + R2L * R2L) / (R2L * R2L);
  // This would reflect q_{maq} to be |q| from page 10/456 of https://github.com/l-korous/doctoral-thesis/blob/master/_reference/Modeling%20of%20H%CE%B1%20Eruptive%20Events%20Observed%20at%20the%20Solar.pdf

  // Sign of winding: corresponds to sign of I_O in TD paper
  iSgn = (td_parameters.N_t >= 0) ? 1.0 : -1.0;

  // "Helicity" factor inside tho loop (used later in B_theta_internal calcs)
  H = 2.0 * (td_parameters.N_t * td_parameters.N_t) / (td_parameters.R * td_parameters.R);
}

/***************************************************************************
Calculate the field according to TD paper (A&A 351, 707, 1999)
Fill the structure with gravity-stratified plasma.
***************************************************************************/
template <EquationsType equationsType, int dim>
void InitialConditionTitovDemoulin<equationsType, dim>::vector_value(const std::vector<Point<dim> > &points, std::vector<std::array<double, Equations<equationsType, dim>::n_components> >& value_list) const
{
  //========== Calculate the vector potential for I_t-generated toroidal field 
  double xx, yy, zz;
  Point<dim> &ca = this->parameters.corner_a;
  Point<dim> &cb = this->parameters.corner_b;

  for (unsigned int pp = 0; pp < points.size(); ++pp)
  {
    const Point<dim>& p = points[pp];

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

    double rho_0 = 0.5 * (1.0 - td_parameters.Tc2Tp) * tanh(r_min - 1.0) + 0.5 * (1 + td_parameters.Tc2Tp);

    if (r_min > 1.0) { // external region

      value_list[pp][0] = rho_0 * exp(-zz * invL_G);         // mass density outside

      B_loc.sadd(1.0, iSgn * td_parameters.R / r_maj, theta0);

      pressure = td_parameters.beta * exp(-zz * invL_G);
    }
    else { // inside the torus

      value_list[pp][0] = rho_0 * exp(-zz * td_parameters.Tc2Tp * invL_G);   // mass density in the loop

      B_loc.sadd(1.0, (iSgn * (sqrt(1.0 + H * (1.0 - r_min * r_min)) + td_parameters.R / r_maj - 1.0)), theta0);

      pressure = td_parameters.beta * exp(-zz * invL_G);
    }

    value_list[pp][1] = 0.0;
    value_list[pp][2] = 0.0;                          // momentum density
    value_list[pp][3] = 0.0;

    // energy density
    value_list[pp][4] = (pressure / (this->parameters.gas_gamma - 1.0)) + 0.5 * (B_loc[0] * B_loc[0] + B_loc[1] * B_loc[1] + B_loc[2] * B_loc[2]);

    value_list[pp][5] = B_loc[0];
    value_list[pp][6] = B_loc[1];                  // magnetic field
    value_list[pp][7] = B_loc[2];
  }
}

template class InitialConditionTitovDemoulin<EquationsTypeMhd, 3>;
