#ifndef _EQUATIONS_MHD_H
#define _EQUATIONS_MHD_H

#include "equations.h"
#include "parameters.h"

template <int dim>
class Equations<EquationsTypeMhd, dim>
{
public:
  // Self-explanatory
  static const unsigned int n_components = 2 * dim + 2;
  
  typedef dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> InputVector;
  
  // Equations constructor takes parameters as an attribute - to set up e.g. gas Gamma value
  Equations(Parameters<dim>& parameters);

  static double compute_kinetic_energy(const InputVector &W);

  static double compute_magnetic_energy(const InputVector &W);

  // Compute pressure, and use kinetic energy and magnetic energy from the state vector.
  double compute_pressure(const InputVector &W) const;

  static double compute_kinetic_energy(const std_cxx11::array<double, n_components> &W);

  static double compute_magnetic_energy(const std_cxx11::array<double, n_components> &W);

  // Compute pressure, and use kinetic energy and magnetic energy from the state vector.
  double compute_pressure(const std_cxx11::array<double, n_components> &W) const;

  // Compute pressure, and use kinetic energy and magnetic energy from the state vector.
  double compute_total_pressure(const InputVector &W) const;

  // Compute pressure, and use the passed values of kinetic energy and magnetic energy.
  double compute_pressure(const InputVector &W, const double& Uk, const double& Um) const;

  double compute_magnetic_field_divergence(const std::vector<Tensor<1, dim> > &W) const;

  void Q(std_cxx11::array<double, n_components> &result, const InputVector &W, const Tensor<1, 3> &normal) const;
  void Q_inv(std_cxx11::array<double, n_components> &result, std_cxx11::array<double, n_components> &F, const Tensor<1, dim> &normal) const;

  // Compute the matrix of MHD fluxes.
  void compute_flux_matrix(const InputVector &W, std_cxx11::array <std_cxx11::array <double, dim>, n_components > &flux) const;

  // Compute the values for the numerical flux
  void numerical_normal_flux(const Tensor<1, dim> &normal, const InputVector &Wplus, const InputVector &Wminus, std_cxx11::array<double, n_components> &normal_flux);

  void store_max_signal_speed(double val);

  void store_max_signal_speed(typename dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double>, false, 1>::value_type val);

  // Passed as a constructor parameter
  Parameters<dim>& parameters;

  // The rest is for the output.  
  class Postprocessor : public DataPostprocessor<dim>
  {
  public:
    Postprocessor(Equations<EquationsTypeMhd, dim>& equations);
    
#if DEAL_II_VERSION_MAJOR > 8 || (DEAL_II_VERSION_MAJOR == 8 && DEAL_II_VERSION_MINOR > 5) || (DEAL_II_VERSION_MAJOR == 8 && DEAL_II_VERSION_MINOR == 5 && DEAL_II_VERSION_SUBMINOR > 0)
    virtual void evaluate_vector_field(const ::DataPostprocessorInputs::Vector<dim> &inputs,
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

  // For CFL.
  double max_signal_speed;

  static std::vector<DataComponentInterpretation::DataComponentInterpretation> component_interpretation();
};

#endif
