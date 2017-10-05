#ifndef _INITIAL_CONDITION_H
#define _INITIAL_CONDITION_H

#include "util.h"
#include "parameters.h"
#include "equations.h"

// Initial condition
template <EquationsType equationsType, int dim>
class InitialCondition : public Function<dim>
{
public:
  // InitialCondition constructor takes parameters as an attribute - to set up whatever necessary and use it in the initial values definition
  InitialCondition(Parameters<dim>&);

  // To be overwritten.
  virtual void vector_value(const std::vector<Point<dim> >& , std::vector<Vector<double> >&) const;
  Parameters<dim>& getParams() const;

private:
  Parameters<dim>& parameters;
};

template <EquationsType equationsType, int dim>
class SimpleICEuler : public InitialCondition<equationsType, dim>
{
public:
  SimpleICEuler(Parameters<dim>&);
  void vector_value(const std::vector<Point<dim> >&, std::vector<Vector<double> >&) const;
};

template <EquationsType equationsType, int dim>
class SimpleICMHD : public InitialCondition<equationsType, dim>
{
public:
  SimpleICMHD(Parameters<dim>&);
  void vector_value(const std::vector<Point<dim> >&, std::vector<Vector<double> >&) const;
};

template <EquationsType equationsType, int dim>
class MHDBlastIC : public InitialCondition<equationsType, dim>
{
public:
  MHDBlastIC(Parameters<dim>&);
  void vector_value(const std::vector<Point<dim> >&, std::vector<Vector<double> >&) const;
};


template <EquationsType equationsType, int dim>
class TitovDemoulinIC : public InitialCondition<equationsType,dim>
{
public:
  TitovDemoulinIC(Parameters<dim>&);
  void vector_value(const std::vector<Point<dim> >&, std::vector<Vector<double> >&) const;
  void point_value(const Point<dim>& , Vector<double>&) const;
  
private:
  double beta;
  double Lg;
  double invLg;
  double N_t;
  double R_t;
  double d2R_t;
  double L2R_t;
  double q_mag;
  double iSgn;
  double heliFactor;
  double Tc2Tp;
  double t_rho;
  double densGrad;
};

#endif
