#include "boundaryConditions.h"
#include "equationsMhd.h"

template<>
BoundaryConditions<EquationsTypeMhd, 3>::BoundaryConditions()
{
  for (unsigned int i = 0; i < max_n_boundaries; ++i)
    for (unsigned int di = 0; di < Equations<EquationsTypeMhd, 3>::n_components; ++di)
      kind[i][di] = Equations<EquationsTypeMhd, 3>::outflow_boundary;
  for (unsigned int di = 0; di < Equations<EquationsTypeMhd, 3>::n_components; ++di)
    kind[0][di] = Equations<EquationsTypeMhd, 3>::inflow_boundary;
};

template <EquationsType equationsType, int dim>
void BoundaryConditions<equationsType, dim>::bc_vector_value(int boundary_id, const std::vector<Point<dim> > &points, std::vector<Vector<double> > & result) const
{
  for (unsigned int i = 0; i < points.size(); ++i)
  {
    result[i][0] = 1.;
    result[i][1] = 0.;
    result[i][2] = 0.;
    result[i][3] = 0.;
    result[i][4] = 0.;
    result[i][5] = 0.;
    result[i][6] = 0.;
    result[i][7] = (points[i][0] > 0.5 ? 10. : 11.) / (1.4 - 1.0) + 0.5;
  }
}

template class BoundaryConditions<EquationsTypeMhd, 3>;