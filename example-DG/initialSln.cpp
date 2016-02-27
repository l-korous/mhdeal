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
  //return p[0] * RHO_EXT;
  return 0.;
}

d InitialSlnMomentumY::value(const Point<DIM> &p)
{
  return 0.;
  //return V2_EXT * RHO_EXT;
}

d InitialSlnMomentumZ::value(const Point<DIM> &p)
{
  return 0.;
  //return V3_EXT * RHO_EXT;
}

d InitialSlnEnergy::value(const Point<DIM> &p)
{
  if(std::sqrt(p(0) * p(0) + p(1) * p(1)) <= .1)
  {
    return 10. / (GAMMA - 1.0) + B1_EXT*B1_EXT+B2_EXT*B2_EXT+B3_EXT*B3_EXT;
  }
  else
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
