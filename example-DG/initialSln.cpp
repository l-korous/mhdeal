#include <iostream>
#include "initialSln.h"

using namespace dealii;
using namespace std;

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
  d x = p[0];
  if(x < 1.5)
    return 0.;
  else if (x < 2.5)
    return x-1.5;
  else if (x < 3.5)
    return 3.5-x;
  else
    return 0.;

  return 1;


  //return p[0] * RHO_INIT;
  //return V1_INIT;
}

d InitialSlnMomentumY::value(const Point<DIM> &p)
{
  return 0;
  //return V2_INIT * RHO_INIT;
}

d InitialSlnMomentumZ::value(const Point<DIM> &p)
{
  return 0;
  //return V3_INIT * RHO_INIT;
}
