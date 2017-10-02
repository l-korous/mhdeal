#include "util.h"
#include "equationsMhd.h"

#define NEGLIGIBLE 1e-12

template <int dim>
Equations<EquationsTypeMhd, dim>::Equations(Parameters<dim>& parameters) : parameters(parameters)
{
}

template <int dim>
std::vector<std::string> Equations<EquationsTypeMhd, dim>::component_names()
{
  return{ "density", "momentum", "momentum", "momentum", "magnetic_field", "magnetic_field", "magnetic_field", "energy" };
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Equations<EquationsTypeMhd, dim>::component_interpretation()
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation;
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
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
    divergence += W[d + 4][d];

  return divergence;
}

template <>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, 3>::compute_magnetic_energy(const InputVector &W)
{
  return 0.5 * (W[4] * W[4] + W[5] * W[5] + W[6] * W[6]);
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeMhd, dim>::compute_pressure(const InputVector &W) const
{
  return (this->parameters.gas_gamma - 1.0) * (W[7] - compute_kinetic_energy(W) - compute_magnetic_energy(W));
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
  return (this->parameters.gas_gamma - 1.0) * (W[7] - Uk - Um);
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::compute_flux_vector(const int derivative, const InputVector &W, std_cxx11::array<typename InputVector::value_type, n_components > &flux) const
{
  const typename InputVector::value_type total_pressure = compute_total_pressure(W);
  typename InputVector::value_type energy_item[dim] = { W[1] * W[4], W[2] * W[5], W[3] * W[6] };

  switch (derivative)
  {
  case 0:
    flux[0] = W[1];
    flux[1] = (W[1] * W[1] / W[0]) - W[4] * W[4] + total_pressure;
    flux[2] = W[1] * W[2] / W[0] - W[4] * W[5];
    flux[3] = W[1] * W[3] / W[0] - W[4] * W[6];
    flux[4] = (W[7] + total_pressure) * W[1] / W[0] - W[4] * energy_item[0];
    flux[5] = 0.0;
    flux[6] = W[1] / W[0] * W[5] - W[2] / W[0] * W[4];
    flux[7] = W[1] / W[0] * W[6] - W[3] / W[0] * W[4];
    break;
  case 1:
    flux[0] = W[2];
    flux[1] = W[2] * W[1] / W[0] - W[5] * W[4];
    flux[2] = W[2] * W[2] / W[0] - W[5] * W[5] + total_pressure;
    flux[3] = W[2] * W[3] / W[0] - W[5] * W[6];
    flux[4] = (W[7] + total_pressure) * W[2] / W[0] - W[5] * energy_item[1];
    flux[5] = W[2] / W[0] * W[4] - W[1] / W[0] * W[5];
    flux[6] = 0.0;
    flux[7] = W[2] / W[0] * W[6] - W[3] / W[0] * W[5];
    break;
  case 2:
    flux[0] = W[3];
    flux[1] = W[3] * W[1] / W[0] - W[6] * W[4];
    flux[2] = W[3] * W[2] / W[0] - W[6] * W[5];
    flux[3] = W[3] * W[3] / W[0] - W[6] * W[6] + total_pressure;
    flux[4] = (W[7] + total_pressure) * W[3] / W[0] - W[6] * energy_item[2];
    flux[5] = W[3] / W[0] * W[4] - W[1] / W[0] * W[6];
    flux[6] = W[3] / W[0] * W[5] - W[2] / W[0] * W[6];
    flux[7] = 0.0;
  }
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::compute_flux_matrix(const InputVector &W, std_cxx11::array <std_cxx11::array <typename InputVector::value_type, dim>, n_components > &flux) const
{
  const typename InputVector::value_type total_pressure = compute_total_pressure(W);
  typename InputVector::value_type energy_item[dim] = { W[1] * W[4], W[2] * W[5], W[3] * W[6] };

  flux[0][0] = W[1];
  flux[1][0] = (W[1] * W[1] / W[0]) - W[4] * W[4] + total_pressure;
  flux[2][0] = W[1] * W[2] / W[0] - W[4] * W[5];
  flux[3][0] = W[1] * W[3] / W[0] - W[4] * W[6];
  flux[4][0] = (W[7] + total_pressure) * W[1] / W[0] - W[4] * energy_item[0];
  flux[5][0] = 0.0;
  flux[6][0] = W[1] / W[0] * W[5] - W[2] / W[0] * W[4];
  flux[7][0] = W[1] / W[0] * W[6] - W[3] / W[0] * W[4];

  flux[0][1] = W[2];
  flux[1][1] = W[2] * W[1] / W[0] - W[5] * W[4];
  flux[2][1] = W[2] * W[2] / W[0] - W[5] * W[5] + total_pressure;
  flux[3][1] = W[2] * W[3] / W[0] - W[5] * W[6];
  flux[4][1] = (W[7] + total_pressure) * W[2] / W[0] - W[5] * energy_item[1];
  flux[5][1] = W[2] / W[0] * W[4] - W[1] / W[0] * W[5];
  flux[6][1] = 0.0;
  flux[7][1] = W[2] / W[0] * W[6] - W[3] / W[0] * W[5];

  flux[0][2] = W[3];
  flux[1][2] = W[3] * W[1] / W[0] - W[6] * W[4];
  flux[2][2] = W[3] * W[2] / W[0] - W[6] * W[5];
  flux[3][2] = W[3] * W[3] / W[0] - W[6] * W[6] + total_pressure;
  flux[4][2] = (W[7] + total_pressure) * W[3] / W[0] - W[6] * energy_item[2];
  flux[5][2] = W[3] / W[0] * W[4] - W[1] / W[0] * W[6];
  flux[6][2] = W[3] / W[0] * W[5] - W[2] / W[0] * W[6];
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

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeMhd, dim>::numerical_normal_flux(const Tensor<1, dim> &normal, const InputVector &Wplus, const InputVector &Wminus,
  std_cxx11::array<typename InputVector::value_type, n_components> &normal_flux)
{
  // If we use Lax Friedrich's flux.
  if (this->parameters.num_flux_type == Parameters<dim>::lax_friedrich)
  {
    std_cxx11::array<std_cxx11::array <typename InputVector::value_type, dim>, n_components > iflux, oflux;

    compute_flux_matrix(Wplus, iflux);
    compute_flux_matrix(Wminus, oflux);

    for (unsigned int di = 0; di < n_components; ++di)
    {
      normal_flux[di] = 0;
      for (unsigned int d = 0; d < dim; ++d)
        normal_flux[di] += 0.5*(iflux[di][d] + oflux[di][d]) * normal[d];

      normal_flux[di] += 0.5*this->parameters.lax_friedrich_stabilization_value*(Wplus[di] - Wminus[di]);
    }

    return;
  }

  // If we use HLLD
  if (this->parameters.num_flux_type == Parameters<dim>::hlld)
  {
    typename InputVector::value_type rho_L, rho_R, VelN_L, VelN_R, MagN_L, MagN_R;
    int dir_abs = (std::abs(normal[0] - 1.) < NEGLIGIBLE ? 0 : (std::abs(normal[1] - 1.) < NEGLIGIBLE ? 1 : 2));
    int dir_sign = normal[dir_abs] > 0. ? 1. : -1.;
    rho_L = Wplus[0], rho_R = Wminus[0];
    VelN_L = Wplus[1 + dir_abs] / rho_L, VelN_R = Wminus[1 + dir_abs] / rho_R;
    MagN_L = Wplus[4 + dir_abs], MagN_R = Wminus[4 + dir_abs];

    typename InputVector::value_type p_L, p_R, p_T_R, p_T_L, S_M, B_x, psT, B2L, B2R, S_L, S_R;
    typename InputVector::value_type aL2, aR2, cfL, cfR;

    B2L = 2. * this->compute_magnetic_energy(Wplus);
    p_L = this->compute_pressure(Wplus);
    p_T_L = p_L + 0.5 * B2L;

    B2R = 2. * this->compute_magnetic_energy(Wminus);
    p_R = this->compute_pressure(Wminus);
    p_T_R = p_R + 0.5 * B2R;

    aL2 = parameters.gas_gamma * p_L / rho_L;
    aR2 = parameters.gas_gamma * p_R / rho_R;

    B2L = B2L / rho_L;
    B2R = B2R / rho_R;
    cfL = std::sqrt(0.5 * (aL2 + B2L + sqrt((aL2 + B2L) * (aL2 + B2L) - 4. * aL2 * MagN_L * MagN_L / rho_L)));
    cfR = std::sqrt(0.5 * (aR2 + B2R + sqrt((aR2 + B2R) * (aR2 + B2R) - 4. * aR2 * MagN_R * MagN_R / rho_R)));

    // Edge states.
    S_L = std::min(VelN_L - cfL, VelN_R - cfR);
    S_R = std::max(VelN_L + cfL, VelN_R + cfR);

    this->store_max_signal_speed(std::max(std::abs(S_L), std::abs(S_R)));

    if (S_L > 0.)
    {
      this->compute_flux_vector(dir_abs, Wplus, normal_flux);
      return;
    }

    if (S_R < 0.)
    {
      this->compute_flux_vector(dir_abs, Wminus, normal_flux);
      return;
    }

    // Now we need to calculate S_M & all other stuff.
    S_M = (((S_R - VelN_R) * (rho_R * VelN_R)) - ((S_L - VelN_L) * (rho_L * VelN_L)) - (p_T_R)+(p_T_L)) / (((S_R - VelN_R) * rho_R) - ((S_L - VelN_L) * rho_L));
    B_x = (S_R * MagN_R - S_L * MagN_L) / (S_R - S_L);

    typename InputVector::value_type scrch1L = S_L - VelN_L, scrch2L = S_L - S_M, scrch3L = S_M - VelN_L;
    typename InputVector::value_type scrch1R = S_R - VelN_R, scrch2R = S_R - S_M, scrch3R = S_M - VelN_R;

    typename InputVector::value_type denom = rho_L * scrch1L * scrch2L - MagN_L * MagN_L;
    typename InputVector::value_type scrch4L = scrch3L / denom, scrch5L = (rho_L * scrch1L * scrch1L - MagN_L * MagN_L) / denom;

    denom = rho_R * scrch1R * scrch2R - MagN_R * MagN_R;
    typename InputVector::value_type scrch4R = scrch3R / denom, scrch5R = (rho_R * scrch1R * scrch1R - MagN_R * MagN_R) / denom;

    psT = (((S_R - VelN_R) * (rho_R * p_T_L)) - ((S_L - VelN_L) * (rho_L * p_T_R)) + ((rho_L * rho_R) * (S_R - VelN_R) * (S_L - VelN_L) * (VelN_R - VelN_L))) / (((S_R - VelN_R) * rho_R) - ((S_L - VelN_L) * rho_L));

    std_cxx11::array<typename InputVector::value_type, n_components> Us_L;
    {
      Us_L[0] = rho_L * (S_L - VelN_L) / (S_L - S_M);
      //////////////////////
      //////////////////////
      //////////////////////
      //////////////////////
    }

    std_cxx11::array<typename InputVector::value_type, n_components> Us_R;
    {
      Us_R[0] = rho_R * (S_R - VelN_R) / (S_R - S_M);
      //////////////////////
      //////////////////////
      //////////////////////
      //////////////////////
    }

    typename InputVector::value_type Ss_L = S_M - (std::abs(B_x) / std::sqrt(Us_L[0]));
    typename InputVector::value_type Ss_R = S_M + (std::abs(B_x) / std::sqrt(Us_R[0]));

    // Fs_L, Fss_L
    if (S_M >= 0.)
    {
      std_cxx11::array<typename InputVector::value_type, n_components> F_L;
      this->compute_flux_vector(dir_abs, Wplus, F_L);
      // THIS is different to RIGHT state.
      if (Ss_L >= 0.)
      {
        for (int k = 0; k < n_components; k++)
          normal_flux[k] = F_L[k] + S_L * (Us_L[k] - Wplus[k]);
        return;
      }
      else
      {
        std_cxx11::array<typename InputVector::value_type, n_components> Uss_L;
        Uss_L[0] = Us_L[0];
        //////////////////////
        //////////////////////
        //////////////////////
        //////////////////////
        for (int k = 0; k < n_components; k++)
          normal_flux[k] = F_L[k] + Ss_L * Uss_L[k] - ((Ss_L - S_L) * Us_L[k]) - (S_L * Wplus[k]);
        return;
      }
    }
    // Fs_R, Fss_R
    else
    {
      std_cxx11::array<typename InputVector::value_type, n_components> F_R;
      this->compute_flux_vector(dir_abs, Wminus, F_R);
      // THIS is different to LEFT state.
      if (Ss_R <= 0.)
      {
        for (int k = 0; k < n_components; k++)
          normal_flux[k] = F_R[k] + S_R * (Us_R[k] - Wminus[k]);
        return;
      }
      else
      {
        std_cxx11::array<typename InputVector::value_type, n_components> Uss_R;
        Uss_R[0] = Us_R[0];

        //////////////////////
        //////////////////////
        //////////////////////
        //////////////////////
        for (int k = 0; k < n_components; k++)
          normal_flux[k] = F_R[k] + Ss_R * Uss_R[k] - ((Ss_R - S_R) * Us_R[k]) - (S_R * Wminus[k]);
        return;
      }
    }
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
  const std::vector<std::vector<Tensor<2, dim> > > &/*dduh*/,
  const std::vector<Point<dim> > &/*normals*/,
  const std::vector<Point<dim> > &/*evaluation_points*/,
  std::vector<Vector<double> > &computed_quantities) const
{
  const unsigned int n_quadrature_points = uh.size();

  Assert(computed_quantities.size() == n_quadrature_points,
    ExcInternalError());

  Assert(uh[0].size() == n_components,
    ExcInternalError());

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
  return update_values;
}

template class Equations<EquationsTypeMhd, 3>;

#include "equationsMhd_inst.cpp"
