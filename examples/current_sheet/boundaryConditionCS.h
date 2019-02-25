#ifndef _BOUNDARY_CONDITION_CS_H
#define _BOUNDARY_CONDITION_CS_H

#include "util.h"
#include "parameters.h"
#include "initialCondition.h"
#include "boundaryCondition.h"
#include "equationsMhd.h"
#include "initialConditionCS.h"
#include "parametersCS.h"

template <int dim>
class BoundaryConditionCSWithVortices : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> values_vector;
  typedef std::array<std::array<double, dim>, Equations<EquationsTypeMhd, dim>::n_components> grad_vector;

  BoundaryConditionCSWithVortices(Parameters<dim>&, CSParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, values_vector &result, 
    const grad_vector &grads, const values_vector &W_plus, double time, typename DoFHandler<dim>::active_cell_iterator&) const;

private:
  CSParameters& cs_parameters;
  double eps;
  double y_1, y_2;
  double r_1_bar(double x, double y) const;
  double r_2_bar(double x, double y) const;
  double omega_1(double x, double y) const;
  double omega_2(double x, double y) const;
  double omega(double time) const;
};

template <int dim>
class BoundaryConditionCSFree : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> values_vector;
  typedef std::array<std::array<double, dim>, Equations<EquationsTypeMhd, dim>::n_components> grad_vector;

  BoundaryConditionCSFree(Parameters<dim>&, CSParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, values_vector &result, 
    const grad_vector &grads, const values_vector &W_plus, double time, typename DoFHandler<dim>::active_cell_iterator&) const;
};

template <int dim>
class BoundaryConditionCSTest : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> values_vector;
  typedef std::array<std::array<double, dim>, Equations<EquationsTypeMhd, dim>::n_components> grad_vector;

  BoundaryConditionCSTest(Parameters<dim>&, CSParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, values_vector &result, 
    const grad_vector &grads, const values_vector &W_plus, double time, typename DoFHandler<dim>::active_cell_iterator&) const;
private:
  CSParameters& cs_parameters;
};

template <int dim>
class BoundaryConditionCSInitialState : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> values_vector;
  typedef std::array<std::array<double, dim>, Equations<EquationsTypeMhd, dim>::n_components> grad_vector;

  BoundaryConditionCSInitialState(Parameters<dim>&, CSParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, values_vector &result,
    const grad_vector &grads, const values_vector &W_plus, double time, typename DoFHandler<dim>::active_cell_iterator&) const;

private:
  CSParameters& cs_parameters;
  InitialConditionCS<EquationsTypeMhd, dim> ic;
  double invL_G;
  double iSgn;
  double d2R;
  double H;
  double L2R;
  double R2L;
  double q_mag;
};
#endif
