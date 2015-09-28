#include "common.h"

d InitialSlnRho::value(const Point<DIM> &p)
{
  return RHO_INIT;
}

d InitialSlnEnergy::value(const Point<DIM> &p)
{
  return E_INIT;
}

d InitialSlnMomentumX::value(const Point<DIM> &p)
{
  return V1_INIT * RHO_INIT;
}

d InitialSlnMomentumY::value(const Point<DIM> &p)
{
  return V2_INIT * RHO_INIT;
}

d InitialSlnMomentumZ::value(const Point<DIM> &p)
{
  return V3_INIT * RHO_INIT;
}