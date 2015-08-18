#include "common.h"

void NumFluxLaxFriedrichs::calculate(vec U_L, vec U_R, Point<DIM> quadPoint, Point<DIM> normal, vec result)
{
  // TODO - this has to be upgraded
  for (ui i = 0; i < COMPONENT_COUNT; i++)
    result[i] = 0.5 * (U_L[i] + U_R[i]);
}

void NumFluxUpwind::calculate(vec U_L, vec U_R, Point<DIM> quadPoint, Point<DIM> normal, vec result)
{
    double norm = std::max<double>(1e-12, std::sqrt(std::pow(quadPoint(0), 2.) + std::pow(quadPoint(1), 2.)));
    d flux = -quadPoint(1) / norm*normal(0) + quadPoint(0) / norm*normal(1);

    for (ui i = 0; i < COMPONENT_COUNT; i++)
        result[i] = flux * (flux >= 0. ? U_L[i] : U_R[i]);
}