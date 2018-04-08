#ifndef _NUM_FLUX_H
#define _NUM_FLUX_H

#include "equationsMhd.h"
#include "parameters.h"

#define n_comp Equations<equationsType, dim>::n_components
#define n_comp_array std::array<double, Equations<equationsType, dim>::n_components>

template <EquationsType equationsType, int dim>
class NumFlux
{
public:
  NumFlux(Parameters<dim>& parameters) : parameters(parameters) {};
  static void Q(n_comp_array &result, const n_comp_array &W, const Tensor<1, dim> &normal);
  static void Q_inv(n_comp_array &result, n_comp_array &F, const Tensor<1, dim> &normal);

  // Compute the values for the numerical flux
  virtual void numerical_normal_flux(const Tensor<1, dim> &normal, const n_comp_array &Wplus,
    const n_comp_array &Wminus, n_comp_array &normal_flux, double& max_speed) const = 0;
protected:
  Parameters<dim>& parameters;
};

template <EquationsType equationsType, int dim>
class NumFluxLaxFriedrich : public NumFlux<equationsType, dim>
{
public:
  NumFluxLaxFriedrich(Parameters<dim>& parameters) : NumFlux<equationsType, dim>(parameters) {};
  void numerical_normal_flux(const Tensor<1, dim> &normal, const n_comp_array &Wplus,
    const n_comp_array &Wminus, n_comp_array &normal_flux, double& max_speed) const;
};

template <EquationsType equationsType, int dim>
class NumFluxHLLD : public NumFlux<equationsType, dim>
{
public:
  NumFluxHLLD(Parameters<dim>& parameters) : NumFlux<equationsType, dim>(parameters) {};
  void numerical_normal_flux(const Tensor<1, dim> &normal, const n_comp_array &Wplus,
    const n_comp_array &Wminus, n_comp_array &normal_flux, double& max_speed) const;
};

#endif
