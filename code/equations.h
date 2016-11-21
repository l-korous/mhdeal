#ifndef _EQUATIONS_H
#define _EQUATIONS_H

enum EquationsType {
  EquationsTypeEuler,
  equationsTypeMhd
};

template <EquationsType equationsType, int dim>
class Equations
{

};

#endif