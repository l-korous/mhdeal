#include "boundaryConditions.h"

using namespace dealii;

d DirichletBoundaryCondition::calculate(ui component, Point<DIM> point)
{
  bool top = (point(0) > 1e-8);

  switch (component)
  {
  case 0:
    if(top)
      return RHO_IN_TOP;
    else
      return RHO_IN_LEFT;
  case 1:
    if (top)
      return RHO_IN_TOP * V1_IN_TOP;
    else
      return RHO_IN_LEFT * V1_IN_LEFT;
  case 2:
    if(top)
      return RHO_IN_TOP * V2_IN_TOP;
    else
      return RHO_IN_LEFT * V2_IN_LEFT;
  case 3:
    if(top)
      return RHO_IN_TOP * V3_IN_TOP;
    else
      return RHO_IN_LEFT * V3_IN_LEFT;
  case 4:
    if(top)
      return E_IN_TOP;
    else
      return E_IN_LEFT;
  }

  abort();
  return 0.;
}
