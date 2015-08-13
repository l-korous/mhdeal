#include "common.h"

void NumFluxLaxFriedrichs::calculate(vec U_L, vec U_R, Point<DIM> normal, vec result)
{
  // TODO - this has to be upgraded
  for (ui i = 0; i < COMPONENT_COUNT; i++)
    result[i] = 0.5 * (U_L[i] + U_R[i]);
}