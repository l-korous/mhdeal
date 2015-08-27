#include "common.h"

void NumFluxLaxFriedrichs::calculate(vec U_L, vec U_R, Point<DIM> quadPoint, Point<DIM> normal, vec& result)
{
  // TODO - this has to be upgraded
  for (ui i = 0; i < COMPONENT_COUNT; i++)
    result[i] = 0.5 * (U_L[i] + U_R[i]);
}

void NumFluxUpwind::calculate(vec U_L, vec U_R, Point<DIM> quadPoint, Point<DIM> normal, vec& result)
{
    d x = quadPoint(0), y = quadPoint(1);
    d nx = normal(0), ny = normal(1);

    for (ui i = 0; i < COMPONENT_COUNT; i++)
        result[i] = FLUX * (FLUX >= 0. ? U_L[i] : U_R[i]);
}