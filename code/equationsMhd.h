#ifndef _EQUATIONS_MHD_H
#define _EQUATIONS_MHD_H

#include "equations.h"
#include "parameters.h"

template <int dim>
class Equations<EquationsTypeMhd, dim>
{
public:
  Equations(Parameters<dim>& parameters);

  static const unsigned int n_components = 2 * dim + 2;
  static const unsigned int first_momentum_component = 1;
  static const unsigned int first_magnetic_flux_component = dim + 1;
  static const unsigned int density_component = 0;
  static const unsigned int energy_component = 2 * dim + 1;

  static std::vector<std::string> component_names();

  static std::vector<DataComponentInterpretation::DataComponentInterpretation> component_interpretation();

  template <typename InputVector>
  typename InputVector::value_type compute_kinetic_energy(const InputVector &W) const;

  template <typename InputVector>
  typename InputVector::value_type compute_magnetic_energy(const InputVector &W) const;

  template <typename InputVector>
  typename InputVector::value_type compute_pressure(const InputVector &W) const;

  template <typename InputVector>
  void compute_flux_matrix(const InputVector &W, std_cxx11::array <std_cxx11::array <typename InputVector::value_type, dim>, n_components > &flux) const;

  template <typename InputVector, typename ValueType>
  void compute_jacobian_addition(double cell_diameter, const InputVector& grad_W, std_cxx11::array <std_cxx11::array <ValueType, dim>, n_components > &jacobian_addition) const;

  template <typename InputVector>
  void numerical_normal_flux(const Tensor<1, dim> &normal, const InputVector &Wplus, const InputVector &Wminus,
    std_cxx11::array<typename InputVector::value_type, n_components> &normal_flux) const;

  template <typename InputVector>
  void compute_forcing_vector(const InputVector &W, std_cxx11::array<typename InputVector::value_type, n_components> &forcing) const;

  template <typename InputVector>
  void Q(std_cxx11::array<typename InputVector::value_type, n_components> &result, const InputVector &W, const Tensor<1, dim> &normal) const;

  template <typename InputVector>
  void Q_inv(std_cxx11::array<typename InputVector::value_type, n_components> &result, std_cxx11::array<typename InputVector::value_type, n_components> &F, const Tensor<1, dim> &normal) const;

  enum BoundaryKind
  {
    inflow_boundary,
    outflow_boundary
  };

  Parameters<dim>& parameters;

  template <typename DataVector>
  void compute_Wminus(const BoundaryKind(&boundary_kind)[n_components], const Tensor<1, dim> &normal_vector, const DataVector &Wplus, const Vector<double> &boundary_values,
    const DataVector &Wminus) const;

  class Postprocessor : public DataPostprocessor<dim>
  {
  public:
    Postprocessor(Equations<EquationsTypeMhd, dim>& equations);

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
  private:
    Equations<EquationsTypeMhd, dim>& equations;
  };
};

#endif