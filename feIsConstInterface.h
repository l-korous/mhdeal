#ifndef _FE_IS_CONSTANT_INTERFACE_H
#define _FE_IS_CONSTANT_INTERFACE_H

#include "util.h"

template <int dim, int spacedim = dim>
class FiniteElementIsConstantInterface
{
public:
  FiniteElementIsConstantInterface() {};

  virtual bool is_constant(const unsigned int i) const = 0;
};

#endif