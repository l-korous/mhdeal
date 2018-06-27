#ifndef _EQUATIONS_MHD_H
#define _EQUATIONS_MHD_H

#include "equations.h"
#include "parameters.h"

template <int dim>
class Equations<EquationsTypeMhd, dim>
{
public:
  static const unsigned int n_components = 2 * dim + 2;

  typedef std::array<double, n_components> values_vector;

  static double compute_kinetic_energy(const values_vector &W);

  static double compute_magnetic_energy(const values_vector &W);

  // Compute pressure, and use kinetic energy and magnetic energy from the state vector.
  static double compute_pressure(const values_vector &W, const Parameters<dim>& parameters);

  // Compute pressure, and use kinetic energy and magnetic energy from the state vector.
  static double compute_total_pressure(const values_vector &W, const Parameters<dim>& parameters);

  // Compute pressure, and use the passed values of kinetic energy and magnetic energy.
  static double compute_pressure(const values_vector &W, const double& Uk, const double& Um, const Parameters<dim>& parameters);

  static double compute_magnetic_field_divergence(const std::vector<Tensor<1, dim> > &W);
  static std::array<double, dim> compute_magnetic_field_curl(const std::vector<Tensor<1, dim> > &W);

  // Compute the matrix of MHD fluxes.
  static void compute_flux_matrix(const values_vector &W, std::array <std::array <double, dim>, n_components > &flux, const Parameters<dim>& parameters);

  // The rest is for the output.  
  class Postprocessor : public DataPostprocessor<dim>
  {
  public:
    Postprocessor(Parameters<dim>& parameters);

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
    Parameters<dim>& parameters;
  };

  static std::vector<std::string> component_names();

  static std::vector<DataComponentInterpretation::DataComponentInterpretation> component_interpretation();
};

#endif
