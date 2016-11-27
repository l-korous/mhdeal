#include "boundaryConditions.h"
#include "equationsEuler.h"
#include "equationsMhd.h"

template<>
BoundaryConditions<EquationsTypeEuler, 2>::BoundaryConditions()
{
  for (unsigned int di = 0; di < Equations<EquationsTypeEuler, 2>::n_components; ++di)
    kind[0][di] = Equations<EquationsTypeEuler, 2>::outflow_boundary;

  for (unsigned int di = 0; di < 2; ++di)
    kind[1][di] = Equations<EquationsTypeEuler, 2>::no_penetration_boundary;

  for (unsigned int di = 2; di < Equations<EquationsTypeEuler, 2>::n_components; ++di)
    kind[1][di] = Equations<EquationsTypeEuler, 2>::outflow_boundary;

  for (unsigned int di = 0; di < Equations<EquationsTypeEuler, 2>::n_components; ++di)
    kind[2][di] = Equations<EquationsTypeEuler, 2>::outflow_boundary;

  for (unsigned int di = 0; di < 2; ++di)
    kind[3][di] = Equations<EquationsTypeEuler, 2>::no_penetration_boundary;

  for (unsigned int di = 2; di < Equations<EquationsTypeEuler, 2>::n_components; ++di)
    kind[3][di] = Equations<EquationsTypeEuler, 2>::outflow_boundary;

  for (unsigned int di = 0; di < 2; ++di)
    kind[4][di] = Equations<EquationsTypeEuler, 2>::no_penetration_boundary;

  for (unsigned int di = 2; di < Equations<EquationsTypeEuler, 2>::n_components; ++di)
    kind[4][di] = Equations<EquationsTypeEuler, 2>::outflow_boundary;
};

template<>
BoundaryConditions<EquationsTypeEuler, 3>::BoundaryConditions()
{
  for (unsigned int di = 0; di < Equations<EquationsTypeEuler, 3>::n_components; ++di)
    kind[0][di] = Equations<EquationsTypeEuler, 3>::outflow_boundary;

  for (unsigned int di = 0; di < 3; ++di)
    kind[1][di] = Equations<EquationsTypeEuler, 3>::no_penetration_boundary;

  for (unsigned int di = 3; di < Equations<EquationsTypeEuler, 3>::n_components; ++di)
    kind[1][di] = Equations<EquationsTypeEuler, 3>::outflow_boundary;

  for (unsigned int di = 0; di < Equations<EquationsTypeEuler, 3>::n_components; ++di)
    kind[2][di] = Equations<EquationsTypeEuler, 3>::outflow_boundary;

  for (unsigned int di = 0; di < 3; ++di)
    kind[3][di] = Equations<EquationsTypeEuler, 3>::no_penetration_boundary;

  for (unsigned int di = 3; di < Equations<EquationsTypeEuler, 3>::n_components; ++di)
    kind[3][di] = Equations<EquationsTypeEuler, 3>::outflow_boundary;

  for (unsigned int di = 0; di < 3; ++di)
    kind[4][di] = Equations<EquationsTypeEuler, 3>::no_penetration_boundary;

  for (unsigned int di = 3; di < Equations<EquationsTypeEuler, 3>::n_components; ++di)
    kind[4][di] = Equations<EquationsTypeEuler, 3>::outflow_boundary;

  for (unsigned int di = 0; di < 3; ++di)
    kind[5][di] = Equations<EquationsTypeEuler, 3>::no_penetration_boundary;

  for (unsigned int di = 3; di < Equations<EquationsTypeEuler, 3>::n_components; ++di)
    kind[5][di] = Equations<EquationsTypeEuler, 3>::outflow_boundary;

  for (unsigned int di = 0; di < 3; ++di)
    kind[6][di] = Equations<EquationsTypeEuler, 3>::no_penetration_boundary;

  for (unsigned int di = 3; di < Equations<EquationsTypeEuler, 3>::n_components; ++di)
    kind[6][di] = Equations<EquationsTypeEuler, 3>::outflow_boundary;
};

template<>
BoundaryConditions<EquationsTypeMhd, 2>::BoundaryConditions()
{
  for (unsigned int i = 0; i < max_n_boundaries; ++i)
    for (unsigned int di = 0; di < Equations<EquationsTypeMhd, 2>::n_components; ++di)
      kind[i][di] = Equations<EquationsTypeMhd, 2>::outflow_boundary;
};

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
      result[i][j] = 0.;
  }
}

template class BoundaryConditions<EquationsTypeEuler, 2>;
template class BoundaryConditions<EquationsTypeEuler, 3>;
template class BoundaryConditions<EquationsTypeMhd, 2>;
template class BoundaryConditions<EquationsTypeMhd, 3>;