#include "util.h"
#include "equations.h"

template <int dim>
std::vector<std::string> Equations<dim>::component_names()
{
  std::vector<std::string> names(dim, "momentum");
  names.push_back("density");
  names.push_back("energy_density");

  return names;
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Equations<dim>::component_interpretation()
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation
    .push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation
    .push_back(DataComponentInterpretation::component_is_scalar);

  return data_component_interpretation;
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<dim>::compute_kinetic_energy(const InputVector &W)
{
  typename InputVector::value_type kinetic_energy = 0;
  for (unsigned int d = 0; d < dim; ++d)
    kinetic_energy += W[first_momentum_component + d] *
    W[first_momentum_component + d];
  kinetic_energy *= 1. / (2 * W[density_component]);

  return kinetic_energy;
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<dim>::compute_pressure(const InputVector &W)
{
  return ((gas_gamma - 1.0) *
    (W[energy_component] - compute_kinetic_energy(W)));
}

template <int dim>
template <typename InputVector>
void Equations<dim>::compute_flux_matrix(const InputVector &W, std_cxx11::array <std_cxx11::array <typename InputVector::value_type, dim>, Equations<dim>::n_components > &flux)
{
  const typename InputVector::value_type pressure = compute_pressure(W);

  for (unsigned int d = 0; d < dim; ++d)
  {
    for (unsigned int e = 0; e < dim; ++e)
      flux[first_momentum_component + d][e] = W[first_momentum_component + d] * W[first_momentum_component + e] / W[density_component];

    flux[first_momentum_component + d][d] += pressure;
  }

  for (unsigned int d = 0; d < dim; ++d)
    flux[density_component][d] = W[first_momentum_component + d];

  for (unsigned int d = 0; d < dim; ++d)
    flux[energy_component][d] = W[first_momentum_component + d] / W[density_component] * (W[energy_component] + pressure);
}

template <int dim>
template <typename InputVector>
void Equations<dim>::numerical_normal_flux(const Tensor<1, dim> &normal, const InputVector &Wplus, const InputVector &Wminus,
  const double alpha, std_cxx11::array<typename InputVector::value_type, n_components> &normal_flux)
{
  std_cxx11::array<std_cxx11::array <typename InputVector::value_type, dim>, Equations<dim>::n_components > iflux, oflux;

  compute_flux_matrix(Wplus, iflux);
  compute_flux_matrix(Wminus, oflux);

  for (unsigned int di = 0; di < n_components; ++di)
  {
    normal_flux[di] = 0;
    for (unsigned int d = 0; d < dim; ++d)
      normal_flux[di] += 0.5*(iflux[di][d] + oflux[di][d]) * normal[d];
   
    normal_flux[di] += 0.5*alpha*(Wplus[di] - Wminus[di]);
  }
}

template <int dim>
template <typename InputVector>
void Equations<dim>::compute_forcing_vector(const InputVector &W, std_cxx11::array<typename InputVector::value_type, n_components> &forcing)
{
  // LK: Tohle u me budou samy nuly
  const double gravity = -1.0;

  for (unsigned int c = 0; c < n_components; ++c)
    switch (c)
    {
    case first_momentum_component + dim - 1:
      forcing[c] = gravity * W[density_component];
      break;
    case energy_component:
      forcing[c] = gravity * W[first_momentum_component + dim - 1];
      break;
    default:
      forcing[c] = 0;
    }
}

template <int dim>
template <typename DataVector>
void Equations<dim>::compute_Wminus(const BoundaryKind(&boundary_kind)[n_components],
  const Tensor<1, dim> &normal_vector,
  const DataVector     &Wplus,
  const Vector<double> &boundary_values,
  const DataVector     &Wminus)
{
  for (unsigned int c = 0; c < n_components; c++)
  {
    switch (boundary_kind[c])
    {
    case inflow_boundary:
    {
      Wminus[c] = boundary_values(c);
      break;
    }

    case outflow_boundary:
    {
      Wminus[c] = Wplus[c];
      break;
    }

    case pressure_boundary:
    {
      const typename DataVector::value_type
        density = (boundary_kind[density_component] ==
          inflow_boundary
          ?
          boundary_values(density_component)
          :
          Wplus[density_component]);

      typename DataVector::value_type kinetic_energy = 0;
      for (unsigned int d = 0; d < dim; ++d)
        if (boundary_kind[d] == inflow_boundary)
          kinetic_energy += boundary_values(d)*boundary_values(d);
        else
          kinetic_energy += Wplus[d] * Wplus[d];
      kinetic_energy *= 1. / 2. / density;

      Wminus[c] = boundary_values(c) / (gas_gamma - 1.0) +
        kinetic_energy;

      break;
    }

    case no_penetration_boundary:
    {
      typename DataVector::value_type vdotn = 0;
      for (unsigned int d = 0; d < dim; d++)
      {
        vdotn += Wplus[d] * normal_vector[d];
      }

      Wminus[c] = Wplus[c] - 2.0*vdotn*normal_vector[c];
      break;
    }

    default:
      Assert(false, ExcNotImplemented());
    }
  }
}

template <int dim>
Equations<dim>::Postprocessor::Postprocessor(const bool do_schlieren_plot) : do_schlieren_plot(do_schlieren_plot)
{}

template <int dim>
void
Equations<dim>::Postprocessor::compute_derived_quantities_vector(
  const std::vector<Vector<double> > &uh,
  const std::vector<std::vector<Tensor<1, dim> > > &duh,
  const std::vector<std::vector<Tensor<2, dim> > > &/*dduh*/,
  const std::vector<Point<dim> > &/*normals*/,
  const std::vector<Point<dim> > &/*evaluation_points*/,
  std::vector<Vector<double> > &computed_quantities) const
{
  const unsigned int n_quadrature_points = uh.size();

  if (do_schlieren_plot == true)
    Assert(duh.size() == n_quadrature_points,
      ExcInternalError());

  Assert(computed_quantities.size() == n_quadrature_points,
    ExcInternalError());

  Assert(uh[0].size() == n_components,
    ExcInternalError());

  if (do_schlieren_plot == true)
    Assert(computed_quantities[0].size() == dim + 2, ExcInternalError())
  else
    Assert(computed_quantities[0].size() == dim + 1, ExcInternalError());

  for (unsigned int q = 0; q < n_quadrature_points; ++q)
  {
    const double density = uh[q](density_component);

    for (unsigned int d = 0; d < dim; ++d)
      computed_quantities[q](d)
      = uh[q](first_momentum_component + d) / density;

    computed_quantities[q](dim) = compute_pressure(uh[q]);

    if (do_schlieren_plot == true)
      computed_quantities[q](dim + 1) = duh[q][density_component] *
      duh[q][density_component];
  }
}

template <int dim>
std::vector<std::string> Equations<dim>::Postprocessor::get_names() const
{
  std::vector<std::string> names;
  for (unsigned int d = 0; d < dim; ++d)
    names.push_back("velocity");
  names.push_back("pressure");

  if (do_schlieren_plot == true)
    names.push_back("schlieren_plot");

  return names;
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Equations<dim>::Postprocessor::get_data_component_interpretation() const
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim,
      DataComponentInterpretation::component_is_part_of_vector);

  interpretation.push_back(DataComponentInterpretation::
    component_is_scalar);

  if (do_schlieren_plot == true)
    interpretation.push_back(DataComponentInterpretation::
      component_is_scalar);

  return interpretation;
}

template <int dim>
UpdateFlags Equations<dim>::Postprocessor::get_needed_update_flags() const
{
  if (do_schlieren_plot == true)
    return update_values | update_gradients;
  else
    return update_values;
}

template class Equations<2>;
template class Equations<3>;

#include "equations_inst.cpp"