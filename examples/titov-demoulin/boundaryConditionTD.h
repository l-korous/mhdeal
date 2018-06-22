#ifndef _BOUNDARY_CONDITION_TD_H
#define _BOUNDARY_CONDITION_TD_H

#include "util.h"
#include "parameters.h"
#include "boundaryCondition.h"
#include "equations.h"
#include "parametersTD.h"

template <int dim>
class BoundaryConditionTDWithVortices : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> InputVector;

  BoundaryConditionTDWithVortices(Parameters<dim>&, TitovDemoulinParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, InputVector &result, const InputVector &W_plus, double time) const;

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
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> InputVector;

  BoundaryConditionTDFree(Parameters<dim>&, TitovDemoulinParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, InputVector &result, const InputVector &W_plus, double time) const;
};

template <int dim>
class BoundaryConditionTDInitialState : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  typedef std::array<double, Equations<EquationsTypeMhd, dim>::n_components> InputVector;

  BoundaryConditionTDInitialState(Parameters<dim>&, TitovDemoulinParameters&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, InputVector &result, const InputVector &W_plus, double time) const;

private:
  TitovDemoulinParameters& td_parameters;
  double invL_G;
  double iSgn;
  double d2R;
  double H;
  double L2R;
  double R2L;
  double q_mag;
};

#endif
