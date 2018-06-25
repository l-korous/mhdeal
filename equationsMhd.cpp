#include "util.h"
#include "equationsMhd.h"

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
inline double Equations<EquationsTypeMhd, dim>::compute_kinetic_energy(const InputVector &W)
{
  return 0.5 * (W[1] * W[1] + W[2] * W[2] + W[3] * W[3]) / W[0];
}

template <int dim>
inline double Equations<EquationsTypeMhd, dim>::compute_magnetic_energy(const InputVector &W)
{
  return 0.5 * (W[5] * W[5] + W[6] * W[6] + W[7] * W[7]);
}

template <int dim>
inline double Equations<EquationsTypeMhd, dim>::compute_pressure(const InputVector &W, const Parameters<dim>& parameters)
{
  return std::max(0., (parameters.gas_gamma - 1.0) * (W[4] - compute_kinetic_energy(W) - compute_magnetic_energy(W)));
}

template <int dim>
inline double Equations<EquationsTypeMhd, dim>::compute_total_pressure(const InputVector &W, const Parameters<dim>& parameters)
{
  const double magnetic_energy = compute_magnetic_energy(W);
  return compute_pressure(W, compute_kinetic_energy(W), magnetic_energy, parameters) + magnetic_energy;
}

template <int dim>
inline double Equations<EquationsTypeMhd, dim>::compute_pressure(const InputVector &W, const double& Uk, const double& Um, const Parameters<dim>& parameters)
{
  return (parameters.gas_gamma - 1.0) * (W[4] - Uk - Um);
}

template <int dim>
double Equations<EquationsTypeMhd, dim>::compute_magnetic_field_divergence(const std::vector<Tensor<1, dim> > &W)
{
  double divergence = 0.;

  for (unsigned int d = 0; d < dim; ++d)
    divergence += W[d + 5][d];

  return divergence;
}

template <int dim>
std::array<double, dim> Equations<EquationsTypeMhd, dim>::compute_magnetic_field_curl(const std::vector<Tensor<1, dim> > &W)
{
  std::array<double, dim> curl_ = { W[7][1] - W[6][2], W[5][2] - W[7][0], W[6][0] - W[5][1] };

  return curl_;
}

template <int dim>
void Equations<EquationsTypeMhd, dim>::compute_flux_matrix(const InputVector &W, std::array <std::array <double, dim>, n_components > &flux, const Parameters<dim>& parameters)
{
  const double mag_energy = compute_magnetic_energy(W);
  const double pressure = compute_pressure(W, parameters);
  const double E = W[4];
  const double total_pressure = pressure + mag_energy;
  const double oneOverRho = 1. / W[0];
  const double UB = (W[1] * W[5] + W[2] * W[6] + W[3] * W[7])* oneOverRho;

  flux[0][0] = W[1];
  flux[1][0] = (W[1] * W[1] * oneOverRho) - W[5] * W[5] + total_pressure;
  flux[2][0] = (W[1] * W[2] * oneOverRho) - W[5] * W[6];
  flux[3][0] = (W[1] * W[3] * oneOverRho) - W[5] * W[7];
  flux[4][0] = (E + total_pressure) * (W[1] * oneOverRho) - (W[5] * UB);
  flux[5][0] = 0.0;
  flux[6][0] = ((W[1] * oneOverRho) * W[6]) - ((W[2] * oneOverRho) * W[5]);
  flux[7][0] = ((W[1] * oneOverRho) * W[7]) - ((W[3] * oneOverRho) * W[5]);

  flux[0][1] = W[2];
  flux[1][1] = (W[2] * W[1] * oneOverRho) - W[6] * W[5];
  flux[2][1] = (W[2] * W[2] * oneOverRho) - W[6] * W[6] + total_pressure;
  flux[3][1] = (W[2] * W[3] * oneOverRho) - W[6] * W[7];
  flux[4][1] = (E + total_pressure) * (W[2] * oneOverRho) - (W[6] * UB);
  flux[5][1] = ((W[2] * oneOverRho) * W[5]) - ((W[1] * oneOverRho) * W[6]);
  flux[6][1] = 0.0;
  flux[7][1] = ((W[2] * oneOverRho) * W[7]) - ((W[3] * oneOverRho) * W[6]);

  flux[0][2] = W[3];
  flux[1][2] = (W[3] * W[1] * oneOverRho) - W[7] * W[5];
  flux[2][2] = (W[3] * W[2] * oneOverRho) - W[7] * W[6];
  flux[3][2] = (W[3] * W[3] * oneOverRho) - W[7] * W[7] + total_pressure;
  flux[4][2] = (E + total_pressure) * (W[3] * oneOverRho) - (W[7] * UB);
  flux[5][2] = ((W[3] * oneOverRho) * W[5]) - ((W[1] * oneOverRho) * W[7]);
  flux[6][2] = ((W[3] * oneOverRho) * W[6]) - ((W[2] * oneOverRho) * W[7]);
  flux[7][2] = 0.0;
}


template <int dim>
Equations<EquationsTypeMhd, dim>::Postprocessor::Postprocessor(Parameters<dim>& parameters) : parameters(parameters)
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
    auto curl_ = equations.compute_magnetic_field_curl(inputs.solution_gradients[q]);
    for (unsigned int d = 0; d < dim; ++d)
      computed_quantities[q](dim + 2 + d) = curl_[d];
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

    std::array<double, n_components> uh_q;
    for (unsigned int d = 0; d < n_components; ++d)
      uh_q[d] = uh[q][d];

    computed_quantities[q](dim) = Equations<EquationsTypeMhd, dim>::compute_pressure(uh_q, parameters);
    computed_quantities[q](dim + 1) = Equations<EquationsTypeMhd, dim>::compute_magnetic_field_divergence(duh[q]);
    auto curl_ = Equations<EquationsTypeMhd, dim>::compute_magnetic_field_curl(duh[q]);
    for (unsigned int d = 0; d < dim; ++d)
      computed_quantities[q](dim + 2 + d) = curl_[d];
  }
}
#endif

template <int dim>
std::vector<std::string> Equations<EquationsTypeMhd, dim>::Postprocessor::get_names() const
{
  return{ "velocity", "velocity", "velocity", "pressure", "divB", "curlB", "curlB", "curlB" };
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Equations<EquationsTypeMhd, dim>::Postprocessor::get_data_component_interpretation() const
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
  interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);

  return interpretation;
}

template <int dim>
UpdateFlags Equations<EquationsTypeMhd, dim>::Postprocessor::get_needed_update_flags() const
{
  return update_values | update_gradients;
}

template class Equations<EquationsTypeMhd, 3>;
