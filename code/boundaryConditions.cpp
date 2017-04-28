#include "boundaryConditions.h"
#include "equationsEuler.h"
#include "equationsMhd.h"

template<>
BoundaryConditions<EquationsTypeMhd, 3>::BoundaryConditions()
{
  for (unsigned int i = 0; i < max_n_boundaries; ++i)
    for (unsigned int di = 0; di < Equations<EquationsTypeMhd, 3>::n_components; ++di)
      kind[i][di] = Equations<EquationsTypeMhd, 3>::outflow_boundary;
};

template <EquationsType equationsType, int dim>
void BoundaryConditions<equationsType, dim>::bc_vector_value(int boundary_id, const std::vector<Point<dim> > &points, std::vector<Vector<double> > & result) const
{
  for (int j = 0; j < Equations<equationsType, dim>::n_components; j++)
  {
    for (int i = 0; i < points.size(); i++)
    {
      result[i][j] = 0.;
    }
  }
}

template class BoundaryConditions<EquationsTypeMhd, 3>;