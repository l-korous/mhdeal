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

d InitialSlnEnergy::value(const Point<DIM> &p)
{
  if (p.norm() < 0.1)
  {
    return 10. / (GAMMA - 1.0) + 1.0;
  }
  else
    return 0.1 / (GAMMA - 1.0) + 1.0;

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
