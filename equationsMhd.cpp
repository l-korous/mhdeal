#include "util.h"
#include "equationsMhd.h"

template <int dim>
Equations<EquationsTypeMhd, dim>::Equations(Parameters<dim>& parameters) : parameters(parameters)
{
}

template <int dim>
std::vector<std::string> Equations<EquationsTypeMhd, dim>::component_names()
{
  return{ "density", "momentum", "momentum", "momentum", "energy", "magnetic_field", "magnetic_field", "magnetic_field" };
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Equations<EquationsTypeMhd, dim>::component_interpretation()
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation;
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  return data_component_interpretation;
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, dim>::compute_kinetic_energy(const InputVector &W)
{
  return 0.5 * (W[1] * W[1] + W[2] * W[2] + W[3] * W[3]) / W[0];
}

template <int dim>
double Equations<EquationsTypeMhd, dim>::compute_magnetic_field_divergence(const std::vector<Tensor<1, dim> > &W) const
{
  double divergence = 0.;

  for (unsigned int d = 0; d < dim; ++d)
    divergence += W[d + 5][d];

  return divergence;
}

template <>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, 3>::compute_magnetic_energy(const InputVector &W)
{
  return 0.5 * (W[5] * W[5] + W[6] * W[6] + W[7] * W[7]);
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, dim>::compute_pressure(const InputVector &W) const
{
  return std::max(0., (this->parameters.gas_gamma - 1.0) * (W[4] - compute_kinetic_energy(W) - compute_magnetic_energy(W)));
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, dim>::compute_total_pressure(const InputVector &W) const
{
  const typename InputVector::value_type magnetic_energy = compute_magnetic_energy(W);
  return compute_pressure(W, compute_kinetic_energy(W), magnetic_energy) + magnetic_energy;
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, dim>::compute_pressure(const InputVector &W, const typename InputVector::value_type& Uk, const typename InputVector::value_type& Um) const
{
  return (this->parameters.gas_gamma - 1.0) * (W[4] - Uk - Um);
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::compute_flux_vector(const int derivative, const InputVector &W, std_cxx11::array<typename InputVector::value_type, n_components > &flux) const
{
  const typename InputVector::value_type mag_energy = compute_magnetic_energy(W);
  const typename InputVector::value_type pressure = compute_pressure(W);
  const typename InputVector::value_type E = W[4] + mag_energy;
  const typename InputVector::value_type total_pressure = pressure + mag_energy;
  const typename InputVector::value_type UB = W[1] * W[5] + W[2] * W[6] + W[3] * W[7];

  switch (derivative)
  {
  case 0:
    flux[0] = W[1];
    flux[1] = (W[1] * W[1] / W[0]) - W[5] * W[5] + total_pressure;
    flux[2] = (W[1] * W[2] / W[0]) - W[5] * W[6];
    flux[3] = (W[1] * W[3] / W[0]) - W[5] * W[7];
    flux[4] = (E + total_pressure) * (W[1] / W[0]) - (W[5] * UB);
    flux[5] = 0.0;
    flux[6] = ((W[1] / W[0]) * W[6]) - ((W[2] / W[0]) * W[5]);
    flux[7] = ((W[1] / W[0]) * W[7]) - ((W[3] / W[0]) * W[5]);
    break;
  case 1:
    flux[0] = W[2];
    flux[1] = (W[2] * W[1] / W[0]) - W[6] * W[5];
    flux[2] = (W[2] * W[2] / W[0]) - W[6] * W[6] + total_pressure;
    flux[3] = (W[2] * W[3] / W[0]) - W[6] * W[7];
    flux[4] = (E + total_pressure) * (W[2] / W[0]) - (W[6] * UB);
    flux[5] = ((W[2] / W[0]) * W[5]) - ((W[1] / W[0]) * W[6]);
    flux[6] = 0.0;
    flux[7] = ((W[2] / W[0]) * W[7]) - ((W[3] / W[0]) * W[6]);
    break;
  case 2:
    flux[0] = W[3];
    flux[1] = (W[3] * W[1] / W[0]) - W[7] * W[5];
    flux[2] = (W[3] * W[2] / W[0]) - W[7] * W[6];
    flux[3] = (W[3] * W[3] / W[0]) - W[7] * W[7] + total_pressure;
    flux[4] = (E + total_pressure) * (W[3] / W[0]) - (W[7] * UB);
    flux[5] = ((W[3] / W[0]) * W[5]) - ((W[1] / W[0]) * W[7]);
    flux[6] = ((W[3] / W[0]) * W[6]) - ((W[2] / W[0]) * W[7]);
    flux[7] = 0.0;
    break;
  }
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::compute_flux_matrix(const InputVector &W, std_cxx11::array <std_cxx11::array <typename InputVector::value_type, dim>, n_components > &flux) const
{
  const typename InputVector::value_type mag_energy = compute_magnetic_energy(W);
  const typename InputVector::value_type pressure = compute_pressure(W);
  const typename InputVector::value_type E = W[4] + mag_energy;
  const typename InputVector::value_type total_pressure = pressure + mag_energy;
  const typename InputVector::value_type UB = W[1] * W[5] + W[2] * W[6] + W[3] * W[7];

  flux[0][0] = W[1];
  flux[1][0] = (W[1] * W[1] / W[0]) - W[5] * W[5] + total_pressure;
  flux[2][0] = (W[1] * W[2] / W[0]) - W[5] * W[6];
  flux[3][0] = (W[1] * W[3] / W[0]) - W[5] * W[7];
  flux[4][0] = (E + total_pressure) * (W[1] / W[0]) - (W[5] * UB);
  flux[5][0] = 0.0;
  flux[6][0] = ((W[1] / W[0]) * W[6]) - ((W[2] / W[0]) * W[5]);
  flux[7][0] = ((W[1] / W[0]) * W[7]) - ((W[3] / W[0]) * W[5]);

  flux[0][1] = W[2];
  flux[1][1] = (W[2] * W[1] / W[0]) - W[6] * W[5];
  flux[2][1] = (W[2] * W[2] / W[0]) - W[6] * W[6] + total_pressure;
  flux[3][1] = (W[2] * W[3] / W[0]) - W[6] * W[7];
  flux[4][1] = (E + total_pressure) * (W[2] / W[0]) - (W[6] * UB);
  flux[5][1] = ((W[2] / W[0]) * W[5]) - ((W[1] / W[0]) * W[6]);
  flux[6][1] = 0.0;
  flux[7][1] = ((W[2] / W[0]) * W[7]) - ((W[3] / W[0]) * W[6]);

  flux[0][2] = W[3];
  flux[1][2] = (W[3] * W[1] / W[0]) - W[7] * W[5];
  flux[2][2] = (W[3] * W[2] / W[0]) - W[7] * W[6];
  flux[3][2] = (W[3] * W[3] / W[0]) - W[7] * W[7] + total_pressure;
  flux[4][2] = (E + total_pressure) * (W[3] / W[0]) - (W[7] * UB);
  flux[5][2] = ((W[3] / W[0]) * W[5]) - ((W[1] / W[0]) * W[7]);
  flux[6][2] = ((W[3] / W[0]) * W[6]) - ((W[2] / W[0]) * W[7]);
  flux[7][2] = 0.0;
}

template <int dim>
template <typename InputVector, typename ValueType>
void Equations<EquationsTypeMhd, dim>::compute_jacobian_addition(double cell_diameter, const InputVector& grad_W, std_cxx11::array <std_cxx11::array <ValueType, dim>, n_components > &jacobian_addition) const
{
  for (unsigned int i = 0; i < n_components; ++i)
    for (unsigned int d = 0; d < dim; ++d)
      jacobian_addition[i][d] = 0.;
}

template <int dim>
void Equations<EquationsTypeMhd, dim>::store_max_signal_speed(typename dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double>, false, 1>::value_type val)
{
  this->max_signal_speed = val.val();
}

template <int dim>
void Equations<EquationsTypeMhd, dim>::store_max_signal_speed(double val)
{
  this->max_signal_speed = val;
}

template <>
template <typename InputVector>
void Equations<EquationsTypeMhd, 3>::Q(std_cxx11::array<typename InputVector::value_type, n_components> &result, const InputVector &W, const Tensor<1, 3> &normal) const
{
  int dir_abs = (std::abs(std::abs(normal[0]) - 1.) < NEGLIGIBLE ? 0 : (std::abs(std::abs(normal[1]) - 1.) < NEGLIGIBLE ? 1 : 2));
  double dir_sign = normal[dir_abs] > 0. ? 1. : -1.;
  result[0] = W[0];
  result[4] = W[4];
  result[1] = W[1 + dir_abs];
  result[5] = W[5 + dir_abs];
  switch (dir_abs) {
  case 0:
    result[2] = W[2];
    result[3] = W[3];
    result[6] = W[6];
    result[7] = W[7];
    break;
  case 1:
    result[2] = W[1];
    result[3] = W[3];
    result[6] = W[5];
    result[7] = W[7];
    break;
  case 2:
    result[2] = W[1];
    result[3] = W[2];
    result[6] = W[5];
    result[7] = W[6];
    break;
  }

  result[1 + dir_abs] *= dir_sign;
  result[5 + dir_abs] *= dir_sign;
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::Q_inv(std_cxx11::array<typename InputVector::value_type, n_components> &result, std_cxx11::array<typename InputVector::value_type, n_components> &W, const Tensor<1, dim> &normal) const
{
  int dir_abs = (std::abs(std::abs(normal[0]) - 1.) < NEGLIGIBLE ? 0 : (std::abs(std::abs(normal[1]) - 1.) < NEGLIGIBLE ? 1 : 2));
  double dir_sign = normal[dir_abs] > 0. ? 1. : -1.;
  result[0] = W[0];
  result[4] = W[4];
  result[1] = W[1 + dir_abs];
  result[5] = W[5 + dir_abs];
  switch (dir_abs) {
  case 0:
    result[2] = W[2];
    result[3] = W[3];
    result[6] = W[6];
    result[7] = W[7];
    break;
  case 1:
    result[1] = W[2];
    result[3] = W[3];
    result[5] = W[6];
    result[7] = W[7];
    break;
  case 2:
    result[1] = W[2];
    result[2] = W[3];
    result[5] = W[6];
    result[6] = W[7];
    break;
  }

  result[1 + dir_abs] *= dir_sign;
  result[5 + dir_abs] *= dir_sign;
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::numerical_normal_flux(const Tensor<1, dim> &normal, const InputVector &Wplus_, const InputVector &Wminus_,
  std_cxx11::array<typename InputVector::value_type, n_components> &normal_flux)
{
  typename InputVector::value_type rho_L, rho_R, VelN_L, VelN_R, MagN_L, MagN_R;
  int dir_abs = (std::abs(std::abs(normal[0]) - 1.) < NEGLIGIBLE ? 0 : (std::abs(std::abs(normal[1]) - 1.) < NEGLIGIBLE ? 1 : 2));
  double dir_sign = normal[dir_abs] > 0. ? 1. : -1.;

  std_cxx11::array <typename InputVector::value_type, n_components> Wplus, Wminus;
  for (unsigned int di = 0; di < 8; ++di)
  {
    Wplus[di] = Wplus_[di];
    Wminus[di] = Wminus_[di];
  }
  Wplus[1 + dir_abs] *= dir_sign;
  Wplus[5 + dir_abs] *= dir_sign;
  Wminus[1 + dir_abs] *= dir_sign;
  Wminus[5 + dir_abs] *= dir_sign;

  // First, calculation of speeds for CFL.
  rho_L = Wplus[0], rho_R = Wminus[0];
  VelN_L = Wplus[1 + dir_abs] / rho_L, VelN_R = Wminus[1 + dir_abs] / rho_R;
  MagN_L = Wplus[5 + dir_abs], MagN_R = Wminus[5 + dir_abs];

  typename InputVector::value_type p_L, p_R, p_T_R, p_T_L, S_M, B_x, psT, BL, BR, B2L, B2R, S_L, S_R;
  typename InputVector::value_type aL2, aR2, cfL, cfR;

  BL = this->compute_magnetic_energy(Wplus);
  p_L = this->compute_pressure(Wplus);
  p_T_L = p_L + BL;

  BR = 2. * this->compute_magnetic_energy(Wminus);
  p_R = this->compute_pressure(Wminus);
  p_T_R = p_R + BR;

  aL2 = parameters.gas_gamma * p_L / rho_L;
  aR2 = parameters.gas_gamma * p_R / rho_R;

  B2L = 2. * BL / rho_L;
  B2R = 2. * BR / rho_R;
  cfL = std::sqrt(0.5 * (aL2 + B2L + sqrt((aL2 + B2L) * (aL2 + B2L) - (4. * aL2 * MagN_L * MagN_L / rho_L))));
  cfR = std::sqrt(0.5 * (aR2 + B2R + sqrt((aR2 + B2R) * (aR2 + B2R) - (4. * aR2 * MagN_R * MagN_R / rho_R))));

  // Edge states.
  S_L = std::min(VelN_L - cfL, VelN_R - cfR);
  S_R = std::max(VelN_L + cfL, VelN_R + cfR);

  this->store_max_signal_speed(std::max(std::abs(S_L), std::abs(S_R)));

  // If we use Lax Friedrich's flux.
  std_cxx11::array<typename InputVector::value_type, n_components> normal_flux_lf;
  {
    std_cxx11::array<std_cxx11::array <typename InputVector::value_type, dim>, n_components > iflux, oflux;

    compute_flux_matrix(Wplus_, iflux);
    compute_flux_matrix(Wminus_, oflux);

    for (unsigned int di = 0; di < n_components; ++di)
    {
      normal_flux[di] = 0;
      for (unsigned int d = 0; d < dim; ++d)
        normal_flux[di] += 0.5 * (iflux[di][d] + oflux[di][d]) * normal[d];

      normal_flux[di] += 0.5 * this->parameters.lax_friedrich_stabilization_value * (Wplus_[di] - Wminus_[di]);
      normal_flux_lf[di] = normal_flux[di];
    }

    if (this->parameters.num_flux_type == Parameters<dim>::lax_friedrich)
      return;
  }

  // If we use HLLD
  if (this->parameters.num_flux_type == Parameters<dim>::hlld)
  {
    typename InputVector::value_type srdl, srdr, hl[2], hr[2], Uk, Um, E2, E3, Sl, Sr, B, B2, cm;
    typename InputVector::value_type sp[5], sml, smr, ptst, ptstr, vbstl, vbstr, Bsgnl, Bsgnr, invsumd;

    std_cxx11::array<typename InputVector::value_type, n_components> Fl;
    std_cxx11::array<typename InputVector::value_type, n_components> Fr;
    std_cxx11::array<typename InputVector::value_type, n_components> Udl;
    std_cxx11::array<typename InputVector::value_type, n_components> Udr;
    std_cxx11::array<typename InputVector::value_type, n_components> Ul;
    std_cxx11::array<typename InputVector::value_type, n_components> Ur;

    std_cxx11::array<typename InputVector::value_type, n_components> QedWminus;
    std_cxx11::array<typename InputVector::value_type, n_components> QedWplus;

    Q(QedWplus, Wplus, normal);
    Q(QedWminus, Wminus, normal);

    // Simple average of mag. field in direction of normal vector
    B = 0.5*(QedWminus[5] + QedWplus[5]);
    B2 = B*B;

    // Calculate left flux
    hl[0] = 1.0 / QedWplus[0];
    Uk = compute_kinetic_energy(QedWplus);
    Um = compute_magnetic_energy(QedWplus);
    hl[1] = compute_pressure(QedWplus, Uk, Um);
    E2 = hl[0] * (QedWplus[1] * QedWplus[7] - QedWplus[3] * QedWplus[5]);
    E3 = hl[0] * (QedWplus[2] * QedWplus[5] - QedWplus[1] * QedWplus[6]);

    Fl[0] = QedWplus[1];
    Fl[1] = hl[0] * QedWplus[1] * QedWplus[1] - QedWplus[5] * QedWplus[5] + Um + hl[1];
    Fl[2] = hl[0] * QedWplus[1] * QedWplus[2] - QedWplus[5] * QedWplus[6];
    Fl[3] = hl[0] * QedWplus[1] * QedWplus[3] - QedWplus[5] * QedWplus[7];
    Fl[5] = 0.0;
    Fl[6] = -E3;
    Fl[7] = E2;
    Fl[4] = hl[0] * QedWplus[1] * (hl[1] * parameters.gas_gamma / (parameters.gas_gamma - 1.0) + Uk) + (E2*QedWplus[7] - E3*QedWplus[6]);

    // Calculate right flux
    hr[0] = 1.0 / QedWminus[0];
    Uk = compute_kinetic_energy(QedWminus);
    Um = compute_magnetic_energy(QedWminus);
    hr[1] = compute_pressure(QedWminus, Uk, Um);
    E2 = hr[0] * (QedWminus[1] * QedWminus[7] - QedWminus[3] * QedWminus[5]);
    E3 = hr[0] * (QedWminus[2] * QedWminus[5] - QedWminus[1] * QedWminus[6]);

    Fr[0] = QedWminus[1];
    Fr[1] = hr[0] * QedWminus[1] * QedWminus[1] - QedWminus[5] * QedWminus[5] + Um + hr[1];
    Fr[2] = hr[0] * QedWminus[1] * QedWminus[2] - QedWminus[5] * QedWminus[6];
    Fr[3] = hr[0] * QedWminus[1] * QedWminus[3] - QedWminus[5] * QedWminus[7];
    Fr[5] = 0.0;
    Fr[6] = -E3;
    Fr[7] = E2;
    Fr[4] = hr[0] * QedWminus[1] * (hr[1] * parameters.gas_gamma / (parameters.gas_gamma - 1.0) + Uk) + (E2*QedWminus[7] - E3*QedWminus[6]);

    // maximum of fast magnetoacoustic speeds L/R
    sp[0] = std::min(QedWplus[1] * hl[0] - cfL, QedWminus[1] * hr[0] - cfR);
    sp[4] = std::max(QedWplus[1] * hl[0] + cfL, QedWminus[1] * hr[0] + cfR);

    // Upwind flux in the case of supersonic flow
    if (sp[0] >= 0.0)
    {
      // use F_L
      for (int k = 0; k < n_components; k++)
        normal_flux[k] = Fl[k];
      Q_inv<InputVector>(normal_flux, normal_flux, normal);
      return;
    }
    if (sp[4] <= 0.0)
    {
      for (int k = 0; k < n_components; k++)
        normal_flux[k] = Fr[k];
      Q_inv<InputVector>(normal_flux, normal_flux, normal);
      return;
    }

    // Determine Alfven and middle speeds
    Sl = sp[0] - QedWplus[1] * hl[0];
    Sr = sp[4] - QedWminus[1] * hr[0];
    sp[2] = (QedWminus[1] * Sr - QedWplus[1] * Sl - p_T_R + p_T_L) / (QedWminus[0] * Sr - QedWplus[0] * Sl);
    sml = sp[0] - sp[2];
    smr = sp[4] - sp[2];

    // Density
    Ul[0] = QedWplus[0] * Sl / sml;
    Ur[0] = QedWminus[0] * Sr / smr;

    Ul[5] = Udl[5] = QedWplus[5];
    Ur[5] = Udr[5] = QedWminus[5];

    srdl = sqrt(Ul[0]);
    srdr = sqrt(Ur[0]);

    // Sl*
    sp[1] = sp[2] - fabs(Ul[5]) / srdl;
    // Sr*
    sp[3] = sp[2] + fabs(Ur[5]) / srdr;

    ptst = p_T_L + QedWplus[0] * Sl*(Sl - sml);
    ptstr = p_T_R + QedWminus[0] * Sr*(Sr - smr);

    // F*_L
    Ul[1] = Ul[0] * sp[2];

    cfL = QedWplus[0] * Sl*sml - Ul[5] * Ul[5];
    if (fabs(cfL) < 1e-8 * ptst)
    {
      Ul[2] = Ul[0] * QedWplus[2] * hl[0];
      Ul[3] = Ul[0] * QedWplus[3] * hl[0];
      Ul[6] = QedWplus[6];
      Ul[7] = QedWplus[7];
    }
    else
    {
      cfL = 1.0 / cfL;
      cm = Ul[5] * (Sl - sml)*cfL;
      Ul[2] = Ul[0] * (QedWplus[2] * hl[0] - QedWplus[6] * cm);
      Ul[3] = Ul[0] * (QedWplus[3] * hl[0] - QedWplus[7] * cm);
      cm = (QedWplus[0] * Sl*Sl - Ul[5] * Ul[5])*cfL;
      Ul[6] = QedWplus[6] * cm;
      Ul[7] = QedWplus[7] * cm;
    }

    vbstl = (Ul[1] * Ul[5] + Ul[2] * Ul[6] + Ul[3] * Ul[7]) / Ul[0];

    Ul[4] = (Sl*QedWplus[4] - p_T_L * QedWplus[1] * hl[0] + ptst * sp[2] + Ul[5] * ((QedWplus[1] * QedWplus[5] + QedWplus[2] * QedWplus[6] + QedWplus[3] * QedWplus[7])*hl[0] - vbstl)) / sml;

    // F*_R
    Ur[1] = Ur[0] * sp[2];
    cfL = QedWminus[0] * Sr*smr - Ur[5] * Ur[5];
    if (fabs(cfL) < 1e-8*ptstr)
    {
      Ur[2] = Ur[0] * QedWminus[2] * hr[0];
      Ur[3] = Ur[0] * QedWminus[3] * hr[0];
      Ur[6] = QedWminus[6];
      Ur[7] = QedWminus[7];
    }
    else
    {
      cfL = 1.0 / cfL;
      cm = Ur[5] * (Sr - smr)*cfL;
      Ur[2] = Ur[0] * (QedWminus[2] * hr[0] - QedWminus[6] * cm);
      Ur[3] = Ur[0] * (QedWminus[3] * hr[0] - QedWminus[7] * cm);
      cm = (QedWminus[0] * Sr*Sr - Ur[5] * Ur[5])*cfL;
      Ur[6] = QedWminus[6] * cm;
      Ur[7] = QedWminus[7] * cm;
    }

    vbstr = (Ur[1] * Ur[5] + Ur[2] * Ur[6] + Ur[3] * Ur[7]) / Ur[0];

    Ur[4] = (Sr*QedWminus[4] - p_T_R * QedWminus[1] * hr[0] + ptstr * sp[2] + Ur[5] * ((QedWminus[1] * QedWminus[5] + QedWminus[2] * QedWminus[6] + QedWminus[3] * QedWminus[7])*hr[0] - vbstr)) / smr;

    if (sp[1] >= 0.0)
    {
      for (register int j = 0; j < n_components; j++)
        normal_flux[j] = Fl[j] + sp[0] * (Ul[j] - QedWplus[j]);

      Q_inv<InputVector>(normal_flux, normal_flux, normal);
      return;
    }
    if (sp[3] <= 0.0 && sp[2] < 0.0)
    {
      for (register int j = 0; j < n_components; j++)
        normal_flux[j] = Fr[j] + sp[4] * (Ur[j] - QedWminus[j]);
      
      Q_inv<InputVector>(normal_flux, normal_flux, normal);
      return;
    }

    // F**_L and F**_R
    if (B2 < 1e-8*(ptst + ptstr))
    {
      for (register int j = 0; j < n_components; j++)
      {
        Udl[j] = Ul[j];
        Udr[j] = Ur[j];
      }
    }
    else
    {
      invsumd = 1.0 / (srdl + srdr);
      Bsgnl = (Ul[5] > 0.0) ? 1.0 : -1.0;
      Bsgnr = (Ur[5] > 0.0) ? 1.0 : -1.0;

      Udl[0] = Ul[0];
      Udr[0] = Ur[0];

      Udl[1] = Ul[1];
      Udr[1] = Ur[1];

      cm = invsumd*(srdl*Ul[2] / Ul[0] + srdr*Ur[2] / Ur[0]);
      cfL = invsumd*(Ur[6] - Ul[6]);
      Udl[2] = Ul[0] * (cm + Bsgnl*cfL);
      Udr[2] = Ur[0] * (cm + Bsgnr*cfL);

      cm = invsumd*(srdl*Ul[3] / Ul[0] + srdr*Ur[3] / Ur[0]);
      cfL = invsumd*(Ur[7] - Ul[7]);
      Udl[3] = Ul[0] * (cm + Bsgnl*cfL);
      Udr[3] = Ur[0] * (cm + Bsgnr*cfL);

      cm = invsumd*(srdl*Ur[6] + srdr*Ul[6]);
      cfL = invsumd*srdl*srdr*(Ur[2] / Ur[0] - Ul[2] / Ul[0]);
      Udl[6] = cm + Bsgnl*cfL;
      Udr[6] = cm + Bsgnr*cfL;

      cm = invsumd*(srdl*Ur[7] + srdr*Ul[7]);
      cfL = invsumd*srdl*srdr*(Ur[3] / Ur[0] - Ul[3] / Ul[0]);
      Udl[7] = cm + Bsgnl*cfL;
      Udr[7] = cm + Bsgnr*cfL;

      Udl[4] = Ul[4] - srdl*Bsgnl*(vbstl - sp[2] * Ul[5] - (Udl[2] * Udl[6] + Udl[3] * Udl[7]) / Udl[0]);
      Udr[4] = Ur[4] + srdr*Bsgnr*(vbstr - sp[2] * Ur[5] - (Udr[2] * Udr[6] + Udr[3] * Udr[7]) / Udr[0]);
    }

    if (sp[2] >= 0.0)
    {
      cm = sp[1] - sp[0];
      for (register int j = 0; j < n_components; j++)
        normal_flux[j] = Fl[j] + sp[1] * Udl[j] - sp[0] * QedWplus[j] - cm*Ul[j];
      
      Q_inv<InputVector>(normal_flux, normal_flux, normal);
    }
    else
    {
      cm = sp[3] - sp[4];
      for (register int j = 0; j < n_components; j++)
        normal_flux[j] = Fr[j] + sp[3] * Udr[j] - sp[4] * QedWminus[j] - cm*Ur[j];
      
      Q_inv<InputVector>(normal_flux, normal_flux, normal);
    }

    return;
  }
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::compute_forcing_vector(const InputVector &W, std_cxx11::array<typename InputVector::value_type, n_components> &forcing) const
{
  for (unsigned int c = 0; c < n_components; ++c)
    forcing[c] = 0.;
}

template <int dim>
Equations<EquationsTypeMhd, dim>::Postprocessor::Postprocessor(Equations<EquationsTypeMhd, dim>& equations) : equations(equations)
{}

#if DEAL_II_VERSION_MAJOR > 8 || (DEAL_II_VERSION_MAJOR == 8 && DEAL_II_VERSION_MINOR > 5) || (DEAL_II_VERSION_MAJOR == 8 && DEAL_II_VERSION_MINOR == 5 && DEAL_II_VERSION_SUBMINOR > 0)
template <int dim>
void
Equations<EquationsTypeMhd, dim>::Postprocessor::evaluate_vector_field(
  const DataPostprocessorInputs::Vector<dim> &inputs,
  std::vector<Vector<double> > &computed_quantities) const
{
  const unsigned int n_quadrature_points = inputs.solution_values.size();

  for (unsigned int q = 0; q < n_quadrature_points; ++q)
  {
    const double density = inputs.solution_values[q](0);

    for (unsigned int d = 0; d < dim; ++d)
      computed_quantities[q](d) = inputs.solution_values[q](1 + d) / density;

    computed_quantities[q](dim) = equations.compute_pressure(inputs.solution_values[q]);
    computed_quantities[q](dim + 1) = equations.compute_magnetic_field_divergence(inputs.solution_gradients[q]);
  }
}
#else
template <int dim>
void
Equations<EquationsTypeMhd, dim>::Postprocessor::compute_derived_quantities_vector(
  const std::vector<Vector<double> > &uh,
  const std::vector<std::vector<Tensor<1, dim> > > &duh,
  const std::vector<std::vector<Tensor<2, dim> > > &,
  const std::vector<Point<dim> > &,
  const std::vector<Point<dim> > &,
  std::vector<Vector<double> > &computed_quantities) const
{
  const unsigned int n_quadrature_points = uh.size();

  for (unsigned int q = 0; q < n_quadrature_points; ++q)
  {
    const double density = uh[q](0);

    for (unsigned int d = 0; d < dim; ++d)
      computed_quantities[q](d) = uh[q](1 + d) / density;

    computed_quantities[q](dim) = equations.compute_pressure(uh[q]);
    computed_quantities[q](dim + 1) = equations.compute_magnetic_field_divergence(duh[q]);
  }
}
#endif

template <int dim>
std::vector<std::string> Equations<EquationsTypeMhd, dim>::Postprocessor::get_names() const
{
  return{ "velocity", "velocity", "velocity", "pressure", "mag_field_divergence" };
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Equations<EquationsTypeMhd, dim>::Postprocessor::get_data_component_interpretation() const
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
  interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  interpretation.push_back(DataComponentInterpretation::component_is_scalar);

  return interpretation;
}

template <int dim>
UpdateFlags Equations<EquationsTypeMhd, dim>::Postprocessor::get_needed_update_flags() const
{
  return update_values | update_gradients;
}

template class Equations<EquationsTypeMhd, 3>;

#include "equationsMhd_inst.cpp"
