#ifndef _NUM_FLUX_H
#define _NUM_FLUX_H

#include "equationsMhd.h"
#include "parameters.h"

template <EquationsType equationsType, int dim>
class NumFlux
{
public:
  static const int n_components = Equations<equationsType, dim>::n_components;
  typedef std::array<double, n_components> component_vector;
  NumFlux(Parameters<dim>& parameters) : parameters(parameters) {};
  static void Q(component_vector &result, const component_vector &W, const Tensor<1, dim> &normal);
  static void Q_inv(component_vector &result, component_vector &F, const Tensor<1, dim> &normal);

  // Compute the values for the numerical flux
  virtual void numerical_normal_flux(const Tensor<1, dim> &normal, const component_vector &Wplus,
    const component_vector &Wminus, component_vector &normal_flux, double& max_speed) const = 0;
protected:
  Parameters<dim>& parameters;
};

template <EquationsType equationsType, int dim>
class NumFluxLaxFriedrich : public NumFlux<equationsType, dim>
{
public:
  NumFluxLaxFriedrich(Parameters<dim>& parameters) : NumFlux<equationsType, dim>(parameters) {};
  void numerical_normal_flux(const Tensor<1, dim> &normal, const component_vector &Wplus,
    const component_vector &Wminus, component_vector &normal_flux, double& max_speed) const;
};

template <EquationsType equationsType, int dim>
class NumFluxHLLD : public NumFlux<equationsType, dim>
{
public:
  NumFluxHLLD(Parameters<dim>& parameters) : NumFlux<equationsType, dim>(parameters) {};
  void numerical_normal_flux(const Tensor<1, dim> &normal, const component_vector &Wplus,
    const component_vector &Wminus, component_vector &normal_flux, double& max_speed) const;
};

#endif
