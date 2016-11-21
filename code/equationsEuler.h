#ifndef _EQUATIONS_EULER_H
#define _EQUATIONS_EULER_H

#include "equations.h"

template <int dim>
class Equations<EquationsTypeEuler, dim>
{
public:
  static const unsigned int n_components = dim + 2;
  static const unsigned int first_momentum_component = 0;
  static const unsigned int density_component = dim;
  static const unsigned int energy_component = dim + 1;

  static std::vector<std::string> component_names();

  static std::vector<DataComponentInterpretation::DataComponentInterpretation> component_interpretation();

  template <typename InputVector>
  static typename InputVector::value_type compute_kinetic_energy(const InputVector &W);

  template <typename InputVector>
  static typename InputVector::value_type compute_pressure(const InputVector &W);

  template <typename InputVector>
  static void compute_flux_matrix(const InputVector &W, std_cxx11::array <std_cxx11::array <typename InputVector::value_type, dim>, Equations<EquationsTypeEuler, dim>::n_components > &flux);

  template <typename InputVector>
  static void numerical_normal_flux(const Tensor<1, dim> &normal, const InputVector &Wplus, const InputVector &Wminus,
    const double alpha, std_cxx11::array<typename InputVector::value_type, n_components> &normal_flux);

  template <typename InputVector>
  static void compute_forcing_vector(const InputVector &W, std_cxx11::array<typename InputVector::value_type, n_components> &forcing);

  enum BoundaryKind
  {
    inflow_boundary,
    outflow_boundary,
    no_penetration_boundary
  };

  template <typename DataVector>
  static void compute_Wminus(const BoundaryKind(&boundary_kind)[n_components], const Tensor<1, dim> &normal_vector, const DataVector &Wplus, const Vector<double> &boundary_values,
    const DataVector &Wminus);

  class Postprocessor : public DataPostprocessor<dim>
  {
  public:
    Postprocessor();

    virtual void compute_derived_quantities_vector(
      const std::vector<Vector<double> > &uh,
      const std::vector<std::vector<Tensor<1, dim> > > &duh,
      const std::vector<std::vector<Tensor<2, dim> > > &dduh,
      const std::vector<Point<dim> > &normals,
      const std::vector<Point<dim> > &evaluation_points,
      std::vector<Vector<double> > &computed_quantities) const;

    virtual std::vector<std::string> get_names() const;

    virtual std::vector<DataComponentInterpretation::DataComponentInterpretation> get_data_component_interpretation() const;

    virtual UpdateFlags get_needed_update_flags() const;
  };
};

#endif