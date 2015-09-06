#include "common.h"

d InitialSlnRho::value(const Point<DIM> &p)
{
  return RHO_EXT;
}

d InitialSlnEnergy::value(const Point<DIM> &p)
{
  return E_EXT;
}

d InitialSlnMomentumX::value(const Point<DIM> &p)
{
  return V1_EXT * RHO_EXT;
}

d InitialSlnMomentumY::value(const Point<DIM> &p)
{
  return V2_EXT * RHO_EXT;
}

d InitialSlnMomentumZ::value(const Point<DIM> &p)
{
  return V3_EXT * RHO_EXT;
}