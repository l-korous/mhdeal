#ifndef _EQUATIONS_H
#define _EQUATIONS_H

#include "util.h"

enum EquationsType {
  EquationsTypeEuler,
  EquationsTypeMhd
};

// Dummy class for templating.
template <EquationsType equationsType, int dim>
class Equations
{
  static const unsigned int n_components = 0;
  typedef std::array<double, n_components> InputVector;
};

#endif