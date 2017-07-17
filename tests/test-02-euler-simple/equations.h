#ifndef _EQUATIONS_H
#define _EQUATIONS_H

enum EquationsType {
  EquationsTypeEuler,
  EquationsTypeMhd
};

// Dummy class for templating.
template <EquationsType equationsType, int dim>
class Equations
{

};

#endif