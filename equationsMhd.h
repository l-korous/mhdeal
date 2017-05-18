#ifndef _EQUATIONS_MHD_H
#define _EQUATIONS_MHD_H

#include "equations.h"
#include "parameters.h"

template <int dim>
class Equations<EquationsTypeMhd, dim>
{
public:
  // Equations constructor takes parameters as an attribute - to set up e.g. gas Gamma value
  Equations(Parameters<dim>& parameters);

  // Self-explanatory
  static const unsigned int n_components = 2 * dim + 2;
  static const unsigned int first_momentum_component = 1;
  static const unsigned int first_magnetic_flux_component = dim + 1;
  static const unsigned int density_component = 0;
  static const unsigned int energy_component = 2 * dim + 1;

  // Boundary types enumeration
  enum BoundaryKind
  {
    inflow_boundary,
    outflow_boundary
  };

  template <typename InputVector>
  typename InputVector::value_type compute_kinetic_energy(const InputVector &W) const;

  double compute_magnetic_field_divergence(const std::vector<Tensor<1, dim> > &W) const;

  template <typename InputVector>
  typename InputVector::value_type compute_magnetic_energy(const InputVector &W) const;

  // Compute pressure, and use kinetic energy and magnetic energy from the state vector.
  template <typename InputVector>
  typename InputVector::value_type compute_pressure(const InputVector &W) const;

  // Compute pressure, and use the passed values of kinetic energy and magnetic energy.
  template <typename InputVector>
  typename InputVector::value_type compute_pressure(const InputVector &W, const typename InputVector::value_type& Uk, const typename InputVector::value_type& Um) const;

  // Compute the matrix of MHD fluxes.
  template <typename InputVector>
  void compute_flux_matrix(const InputVector &W, std_cxx11::array <std_cxx11::array <typename InputVector::value_type, dim>, n_components > &flux) const;

  // Compute jacobian addition - for grad grad | grad div | grad rot terms
  template <typename InputVector, typename ValueType>
  void compute_jacobian_addition(double cell_diameter, const InputVector& grad_W, std_cxx11::array <std_cxx11::array <ValueType, dim>, n_components > &jacobian_addition) const;

  // Compute the values for the numerical flux
  template <typename InputVector>
  void numerical_normal_flux(const Tensor<1, dim> &normal, const InputVector &Wplus, const InputVector &Wminus,
    std_cxx11::array<typename InputVector::value_type, n_components> &normal_flux) const;

  // Helper used in numerical_normal_flux
  // Rotational matrix taking any normal, and rotating the solution so that the new x-coordinate is in the direction of the normal.
  template <typename InputVector>
  void Q(std_cxx11::array<typename InputVector::value_type, n_components> &result, const InputVector &W, const Tensor<1, dim> &normal) const;

  // Helper used in numerical_normal_flux
  // Inverse to Q()
  template <typename InputVector>
  void Q_inv(std_cxx11::array<typename InputVector::value_type, n_components> &result, std_cxx11::array<typename InputVector::value_type, n_components> &F, const Tensor<1, dim> &normal) const;

  // Helper used in numerical_normal_flux
  // Calculate Wminus (state vector on the other-external side of the boundary)
  // E.g. in free BC, this is the same state vector as the provided Wplus, for Dirichlet, these are the prescribed values.
  template <typename InputVector>
  void compute_Wminus(const BoundaryKind(&boundary_kind)[n_components], const Tensor<1, dim> &normal_vector, const InputVector &Wplus, const Vector<double> &boundary_values,
    const InputVector &Wminus) const;

  // Compute the values for forcing (source, absolute) term
  template <typename InputVector>
  void compute_forcing_vector(const InputVector &W, std_cxx11::array<typename InputVector::value_type, n_components> &forcing) const;

  // Passed as a constructor parameter
  Parameters<dim>& parameters;

  // The rest is for the output.  
  class Postprocessor : public DataPostprocessor<dim>
  {
  public:
    Postprocessor(Equations<EquationsTypeMhd, dim>& equations);
    
#if DEAL_II_VERSION_MAJOR > 8 || (DEAL_II_VERSION_MAJOR == 8 && DEAL_II_VERSION_MINOR > 5) || (DEAL_II_VERSION_MAJOR == 8 && DEAL_II_VERSION_MINOR == 5 && DEAL_II_VERSION_SUBMINOR > 0)
    virtual void evaluate_vector_field(
      const dealii::DataPostprocessor::DataPostprocessorInputs::Vector<dim> &inputs,
      std::vector<Vector<double> > &computed_quantities) const;
#else
    virtual void compute_derived_quantities_vector(
      const std::vector<Vector<double> > &uh,
      const std::vector<std::vector<Tensor<1, dim> > > &duh,
      const std::vector<std::vector<Tensor<2, dim> > > &dduh,
      const std::vector<Point<dim> > &normals,
      const std::vector<Point<dim> > &evaluation_points,
      std::vector<Vector<double> > &computed_quantities) const;
#endif

    virtual std::vector<std::string> get_names() const;

    virtual std::vector<DataComponentInterpretation::DataComponentInterpretation> get_data_component_interpretation() const;

    virtual UpdateFlags get_needed_update_flags() const;
  private:
    Equations<EquationsTypeMhd, dim>& equations;
  };

  static std::vector<std::string> component_names();
  static std::vector<DataComponentInterpretation::DataComponentInterpretation> component_interpretation();
};

#endif
