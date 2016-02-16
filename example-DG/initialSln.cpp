#include <iostream>
#include "initialSln.h"

using namespace dealii;
using namespace std;

d InitialSlnRho::value(const Point<DIM> &p)
{
  return RHO_EXT;
}

d InitialSlnMomentumX::value(const Point<DIM> &p)
{
  return p[0] * RHO_EXT;
  return V1_EXT;
}

d InitialSlnMomentumY::value(const Point<DIM> &p)
{
  return 0;
  //return V2_EXT * RHO_EXT;
}

d InitialSlnMomentumZ::value(const Point<DIM> &p)
{
  return 0;
  //return V3_EXT * RHO_EXT;
}

d InitialSlnEnergy::value(const Point<DIM> &p)
{
  return E_EXT;
}

d InitialSlnB1::value(const Point<DIM> &p)
{
  return B1_EXT;
}

d InitialSlnB2::value(const Point<DIM> &p)
{
  return B2_EXT;
}

d InitialSlnB3::value(const Point<DIM> &p)
{
  return B3_EXT;
}

d InitialSlnJ1::value(const Point<DIM> &p)
{
  return J1_EXT;
}

d InitialSlnJ2::value(const Point<DIM> &p)
{
  return J2_EXT;
}

d InitialSlnJ3::value(const Point<DIM> &p)
{
  return J3_EXT;
}
