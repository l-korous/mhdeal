#include "util.h"
#include "equationsEuler.h"

template <int dim>
Equations<EquationsTypeEuler, dim>::Equations(Parameters<dim>& parameters) : parameters(parameters)
{
}

template <int dim>
std::vector<std::string> Equations<EquationsTypeEuler, dim>::component_names()
{
  std::vector<std::string> names(dim, "momentum");
  names.push_back("density");
  names.push_back("energy_density");

  return names;
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Equations<EquationsTypeEuler, dim>::component_interpretation()
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);

  return data_component_interpretation;
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeEuler, dim>::compute_kinetic_energy(const InputVector &W) const
{
  typename InputVector::value_type kinetic_energy = 0;
  
  for (unsigned int d = 0; d < dim; ++d)
    kinetic_energy += W[first_momentum_component + d] * W[first_momentum_component + d];

  kinetic_energy *= 1. / (2 * W[density_component]);

  return kinetic_energy;
}

template <int dim>
template <typename InputVector>
typename InputVector::value_type Equations<EquationsTypeEuler, dim>::compute_pressure(const InputVector &W) const
{
  return ((this->parameters.gas_gamma - 1.0) * (W[energy_component] - compute_kinetic_energy(W)));
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeEuler, dim>::compute_flux_matrix(const InputVector &W, std_cxx11::array <std_cxx11::array <typename InputVector::value_type, dim>, n_components > &flux) const
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
template <typename InputVector, typename ValueType>
void Equations<EquationsTypeEuler, dim>::compute_jacobian_addition(double cell_diameter, const InputVector& grad_W, std_cxx11::array <std_cxx11::array <ValueType, dim>, n_components > &jacobian_addition) const
{
  for (unsigned int i = 0; i < n_components; ++i)
    for (unsigned int d = 0; d < dim; ++d)
      jacobian_addition[i][d] = std::pow(cell_diameter, 2.0) *  grad_W[i][d];
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeEuler, dim>::numerical_normal_flux(const Tensor<1, dim> &normal, const InputVector &Wplus, const InputVector &Wminus,
  std_cxx11::array<typename InputVector::value_type, n_components> &normal_flux) const
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
}

template <int dim>
template <typename InputVector>
void Equations<EquationsTypeEuler, dim>::compute_forcing_vector(const InputVector &W, std_cxx11::array<typename InputVector::value_type, n_components> &forcing) const
{
  // LK: Tohle u me budou samy nuly
  const double gravity = -1.0;

  for (unsigned int c = 0; c < n_components; ++c)
    switch (c)
    {
    case first_momentum_component + 1:
      forcing[c] = gravity * W[density_component];
      break;
    case energy_component:
      forcing[c] = gravity * W[first_momentum_component + 1];
      break;
    default:
      forcing[c] = 0;
    }
}

template <int dim>
template <typename DataVector>
void Equations<EquationsTypeEuler, dim>::compute_Wminus(const BoundaryKind(&boundary_kind)[n_components],
  const Tensor<1, dim> &normal_vector,
  const DataVector     &Wplus,
  const Vector<double> &boundary_values,
  const DataVector     &Wminus) const
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
Equations<EquationsTypeEuler, dim>::Postprocessor::Postprocessor(Equations<EquationsTypeEuler, dim>& equations) : equations(equations)
{}

template <int dim>
void
Equations<EquationsTypeEuler, dim>::Postprocessor::compute_derived_quantities_vector(
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

    Assert(computed_quantities[0].size() == dim + 1, ExcInternalError());

  for (unsigned int q = 0; q < n_quadrature_points; ++q)
  {
    const double density = uh[q](density_component);

    for (unsigned int d = 0; d < dim; ++d)
      computed_quantities[q](d) = uh[q](first_momentum_component + d) / density;

    computed_quantities[q](dim) = equations.compute_pressure(uh[q]);
  }
}

template <int dim>
std::vector<std::string> Equations<EquationsTypeEuler, dim>::Postprocessor::get_names() const
{
  std::vector<std::string> names;
  for (unsigned int d = 0; d < dim; ++d)
    names.push_back("velocity");
  names.push_back("pressure");

  return names;
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Equations<EquationsTypeEuler, dim>::Postprocessor::get_data_component_interpretation() const
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim,
      DataComponentInterpretation::component_is_part_of_vector);

  interpretation.push_back(DataComponentInterpretation::
    component_is_scalar);

  return interpretation;
}

template <int dim>
UpdateFlags Equations<EquationsTypeEuler, dim>::Postprocessor::get_needed_update_flags() const
{
    return update_values;
}

template class Equations<EquationsTypeEuler, 2>;
template class Equations<EquationsTypeEuler, 3>;

#include "equationsEuler_inst.cpp"