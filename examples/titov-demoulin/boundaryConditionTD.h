#ifndef _BOUNDARY_CONDITION_TD_H
#define _BOUNDARY_CONDITION_TD_H

#include "util.h"
#include "parameters.h"
#include "boundaryCondition.h"
#include "equations.h"
#include "initialConditionTD.h"
#include "parametersTD.h"

template <int dim>
class BoundaryConditionTDWithVortices : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> values_vector;
  typedef std::array<std::array<double, dim>, Equations<EquationsTypeMhd, dim>::n_components> grad_vector;

  BoundaryConditionTDWithVortices(Parameters<dim>&, TitovDemoulinParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, values_vector &result, 
    const grad_vector &grads, const values_vector &W_plus, double time, typename DoFHandler<dim>::active_cell_iterator&) const;

private:
  TitovDemoulinParameters& td_parameters;
  double eps;
  double y_1, y_2;
  double r_1_bar(double x, double y) const;
  double r_2_bar(double x, double y) const;
  double omega_1(double x, double y) const;
  double omega_2(double x, double y) const;
  double omega(double time) const;
};

template <int dim>
class BoundaryConditionTDFree : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> values_vector;
  typedef std::array<std::array<double, dim>, Equations<EquationsTypeMhd, dim>::n_components> grad_vector;

  BoundaryConditionTDFree(Parameters<dim>&, TitovDemoulinParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, values_vector &result, 
    const grad_vector &grads, const values_vector &W_plus, double time, typename DoFHandler<dim>::active_cell_iterator&) const;
};

template <int dim>
class BoundaryConditionTDTest : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> values_vector;
  typedef std::array<std::array<double, dim>, Equations<EquationsTypeMhd, dim>::n_components> grad_vector;

  BoundaryConditionTDTest(Parameters<dim>&, TitovDemoulinParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, values_vector &result, 
    const grad_vector &grads, const values_vector &W_plus, double time, typename DoFHandler<dim>::active_cell_iterator&) const;
};

template <int dim>
class BoundaryConditionTDInitialState : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> values_vector;
  typedef std::array<std::array<double, dim>, Equations<EquationsTypeMhd, dim>::n_components> grad_vector;

  BoundaryConditionTDInitialState(Parameters<dim>&, TitovDemoulinParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, values_vector &result,
    const grad_vector &grads, const values_vector &W_plus, double time, typename DoFHandler<dim>::active_cell_iterator&) const;

private:
  TitovDemoulinParameters& td_parameters;
  InitialConditionTitovDemoulin<EquationsTypeMhd, dim> ic;
  double invL_G;
  double iSgn;
  double d2R;
  double H;
  double L2R;
  double R2L;
  double q_mag;
};

#endif
