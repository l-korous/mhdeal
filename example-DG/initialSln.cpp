#include "common.h"

d InitialSlnRho::value(const Point<DIM> &p)
{
  return RHO_EXT;
}

d InitialSlnEnergy::value(const Point<DIM> &p)
{
  return E_EXT;
}
