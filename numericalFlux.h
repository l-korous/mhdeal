#ifndef _NUM_FLUX_H
#define _NUM_FLUX_H

#include "equationsMhd.h"
#include "parameters.h"

template <EquationsType equationsType, int dim>
class NumFlux
{
public:
  NumFlux(Parameters<dim>& parameters) : parameters(parameters) {};
  static void Q(std::array<double, Equations<equationsType, dim>::n_components> &result, const std::array<double, Equations<equationsType, dim>::n_components> &W, const Tensor<1, dim> &normal);
  static void Q_inv(std::array<double, Equations<equationsType, dim>::n_components> &result, std::array<double, Equations<equationsType, dim>::n_components> &F, const Tensor<1, dim> &normal);

  // Compute the values for the numerical flux
  virtual void numerical_normal_flux(const Tensor<1, dim> &normal, const std::array<double, Equations<equationsType, dim>::n_components> &Wplus,
    const std::array<double, Equations<equationsType, dim>::n_components> &Wminus, std::array<double, Equations<equationsType, dim>::n_components> &normal_flux, double& max_speed) const = 0;
protected:
  Parameters<dim>& parameters;
};

template <EquationsType equationsType, int dim>
class NumFluxLaxFriedrich : public NumFlux<equationsType, dim>
{
public:
  NumFluxLaxFriedrich(Parameters<dim>& parameters) : NumFlux<equationsType, dim>(parameters) {};
  void numerical_normal_flux(const Tensor<1, dim> &normal, const std::array<double, Equations<equationsType, dim>::n_components> &Wplus,
    const std::array<double, Equations<equationsType, dim>::n_components> &Wminus, std::array<double, Equations<equationsType, dim>::n_components> &normal_flux, double& max_speed) const;
};

template <EquationsType equationsType, int dim>
class NumFluxHLLD : public NumFlux<equationsType, dim>
{
public:
  NumFluxHLLD(Parameters<dim>& parameters) : NumFlux<equationsType, dim>(parameters) {};
  void numerical_normal_flux(const Tensor<1, dim> &normal, const std::array<double, Equations<equationsType, dim>::n_components> &Wplus,
    const std::array<double, Equations<equationsType, dim>::n_components> &Wminus, std::array<double, Equations<equationsType, dim>::n_components> &normal_flux, double& max_speed) const;
};

#endif
