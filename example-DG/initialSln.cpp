#include "common.h"

d InitialSln::value(const Point<DIM> &p)
{
  d x = p(0), y = p(1);
  d result = std::sin(1. / (x + .1)) * std::cos(1. / (y + .08));
  return result;
}