#include "feDivFree.h"

#define DEGREE 1

template <int dim, int spacedim>
FE_DG_DivFree<dim, spacedim>::FE_DG_DivFree()
  :
  FiniteElement<dim, spacedim>(FiniteElementData<dim>(get_dpo_vector(DEGREE), dim, DEGREE, FiniteElementData<dim>::L2),
    std::vector<bool>(FiniteElementData<dim>(get_dpo_vector(DEGREE), dim, DEGREE).dofs_per_cell, true),
    std::vector<ComponentMask>({ ComponentMask({1, 0, 0}), ComponentMask({0, 1, 0}), ComponentMask({0, 0, 1}), ComponentMask({ 1, 0, 0 }), ComponentMask({ 1, 0, 0 }), ComponentMask({ 0, 1, 0 }), ComponentMask({ 0, 1, 0 }), ComponentMask({ 0, 0, 1 }), ComponentMask({ 0, 0, 1 }), ComponentMask({1, 1, 1}) }))
  , FiniteElementIsConstantInterface<dim>()
  //std::vector<ComponentMask>(10, ComponentMask({ true, true, true })))
{
  this->system_to_component_table[0] = std::pair<unsigned int, unsigned int>(0, 0);
  this->system_to_component_table[1] = std::pair<unsigned int, unsigned int>(1, 1);
  this->system_to_component_table[2] = std::pair<unsigned int, unsigned int>(2, 2);
  this->system_to_component_table[3] = std::pair<unsigned int, unsigned int>(0, 3);
  this->system_to_component_table[4] = std::pair<unsigned int, unsigned int>(0, 4);
  this->system_to_component_table[5] = std::pair<unsigned int, unsigned int>(1, 5);
  this->system_to_component_table[6] = std::pair<unsigned int, unsigned int>(1, 6);
  this->system_to_component_table[7] = std::pair<unsigned int, unsigned int>(2, 7);
  this->system_to_component_table[8] = std::pair<unsigned int, unsigned int>(2, 8);
  const unsigned int n_dofs = this->dofs_per_cell;
  for (unsigned int ref_case = RefinementCase<dim>::cut_x; ref_case < RefinementCase<dim>::isotropic_refinement + 1; ++ref_case)
  {
    // do nothing, as anisotropic refinement is not implemented so far
    if (dim != 2 && ref_case != RefinementCase<dim>::isotropic_refinement)
      continue;

    const unsigned int nc = GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));
    for (unsigned int i = 0; i < nc; ++i)
    {
      this->prolongation[ref_case - 1][i].reinit(n_dofs, n_dofs);
      // Fill prolongation matrices with embedding operators
      for (unsigned int j = 0; j < n_dofs; ++j)
        this->prolongation[ref_case - 1][i](j, j) = 1.;
    }
  }
}


template <int dim, int spacedim>
FiniteElement<dim, spacedim> * FE_DG_DivFree<dim, spacedim>::clone() const
{
  return new FE_DG_DivFree<dim, spacedim>(*this);
}


template <int dim, int spacedim>
std::vector<unsigned int>
FE_DG_DivFree<dim, spacedim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim + 1, static_cast<unsigned int>(0));
  dpo[dim] = 10;
  return dpo;
}


template <int dim, int spacedim>
UpdateFlags FE_DG_DivFree<dim, spacedim>::requires_update_flags(const UpdateFlags flags) const
{
  UpdateFlags out = flags;

  if (flags & (update_values | update_gradients | update_hessians))
    out |= update_quadrature_points;

  return out;
}


template <int dim, int spacedim>
typename FiniteElement<dim, spacedim>::InternalDataBase *
FE_DG_DivFree<dim, spacedim>::get_data(const UpdateFlags update_flags, const Mapping<dim, spacedim> &, const Quadrature<dim> &, dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &/*output_data*/) const
{
  typename FiniteElement<dim, spacedim>::InternalDataBase *data = new typename FiniteElement<dim, spacedim>::InternalDataBase;
  data->update_each = requires_update_flags(update_flags);
  return data;
}


template <int dim, int spacedim>
double FE_DG_DivFree<dim, spacedim>::shape_value(const unsigned int i, const Point<dim> &p) const
{
  Assert(false, ExcNotImplemented());
  return 0;
}


template <int dim, int spacedim>
double
FE_DG_DivFree<dim, spacedim>::shape_value(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const unsigned int i, const Point<dim> &p) const

{
  // Assert(this->is_primitive(), ExcShapeFunctionNotPrimitive());
  Assert(false, ExcNotImplemented());
  return 0.;
}


template <int dim, int spacedim>
void
FE_DG_DivFree<dim, spacedim>::shape_values(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const std::vector< Point<dim> > &p, std::vector< std::vector<double> >& values) const
{
  // Assert(this->is_primitive(), ExcShapeFunctionNotPrimitive());
  Assert(false, ExcNotImplemented());
}


template <int dim, int spacedim>
double FE_DG_DivFree<dim, spacedim>::shape_value_component(const unsigned int i, const Point<dim> &p, const unsigned int component) const
{
  Assert(false, ExcNotImplemented());
  return 0;
}


template <int dim, int spacedim>
double FE_DG_DivFree<dim, spacedim>::shape_value_component(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const unsigned int i, const Point<dim> &p, const unsigned int component) const
{
  Assert(false, ExcNotImplemented());
  return 0;
}


template <int dim, int spacedim>
inline double
FE_DG_DivFree<dim, spacedim>::shape_value_component(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const unsigned int i, const Point<dim> &p, const unsigned int component, double cell_diameter) const
{
  Assert(i < this->dofs_per_cell, ExcIndexRange(i, 0, this->dofs_per_cell));
  switch (i) {
  case 0:
    switch (component) {
    case 0:
      return 1.;
      break;
    }
    break;
  case 1:
    switch (component) {
    case 1:
      return 1.;
      break;
    }
    break;
  case 2:
    switch (component) {
    case 2:
      return 1.;
      break;
    }
    break;
  case 3:
    switch (component) {
    case 0:
      return p(1);
      break;
    }
    break;
  case 4:
    switch (component) {
    case 0:
      return p(2);
      break;
    }
    break;
  case 5:
    switch (component) {
    case 1:
      return p(0);
      break;
    }
    break;
  case 6:
    switch (component) {
    case 1:
      return p(2);
      break;
    }
    break;
  case 7:
    switch (component) {
    case 2:
      return p(0);
      break;
    }
    break;
  case 8:
    switch (component) {
    case 2:
      return p(1);
      break;
    }
    break;
  case 9:
    switch (component) {
    case 0:
      return 2. * cell_diameter * p(0);
      break;
    case 1:
      return -cell_diameter * p(1);
      break;
    case 2:
      return -cell_diameter * p(2);
      break;
    }
    break;
  }
  return 0.;
}


template <int dim, int spacedim>
Tensor<1, dim> FE_DG_DivFree<dim, spacedim>::shape_grad(const unsigned int i, const Point<dim> &p) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<1, dim>();
}


template <int dim, int spacedim>
Tensor<1, dim>
FE_DG_DivFree<dim, spacedim>::shape_grad(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const unsigned int i, const Point<dim> &p) const
{
  // Assert(this->is_primitive(), ExcShapeFunctionNotPrimitive());
  Assert(false, ExcNotImplemented());
  return Tensor<1, dim>();
}


template <int dim, int spacedim>
Tensor<1, dim>
FE_DG_DivFree<dim, spacedim>::shape_grad_component(const unsigned int i, const Point<dim> &p, const unsigned int component) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<1, dim>();
}


template <int dim, int spacedim>
inline Tensor<1, dim>
FE_DG_DivFree<dim, spacedim>::shape_grad_component(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const unsigned int i, const Point<dim> &p, const unsigned int component) const
{
  // TODO: Use extent_in_direction() & mapping instead of diameter()
  const double h = cell->diameter();
  Tensor<1, dim> grad({ 0., 0., 0. });
  switch (i) {
  case 3:
    switch (component) {
    case 0:
      grad[1] = 1. / h;
      break;
    }
    break;
  case 4:
    switch (component) {
    case 0:
      grad[2] = 1. / h;
      break;
    }
    break;
  case 5:
    switch (component) {
    case 1:
      grad[0] = 1. / h;
      break;
    }
    break;
  case 6:
    switch (component) {
    case 1:
      grad[2] = 1. / h;
      break;
    }
    break;
  case 7:
    switch (component) {
    case 2:
      grad[0] = 1. / h;
      break;
    }
    break;
  case 8:
    switch (component) {
    case 2:
      grad[1] = 1. / h;
      break;
    }
    break;
  case 9:
    switch (component) {
    case 0:
      grad[0] = 2.;
      break;
    case 1:
      grad[1] = -1.;
      break;
    case 2:
      grad[2] = -1.;
      break;
    }
    break;
  }
  return grad;
}


template <int dim, int spacedim>
Tensor<2, dim>
FE_DG_DivFree<dim, spacedim>::shape_grad_grad(const unsigned int i, const Point<dim> &p) const
{
  // Assert(this->is_primitive(), ExcShapeFunctionNotPrimitive());
  Assert(false, ExcNotImplemented());
  return Tensor<2, dim>();
}


template <int dim, int spacedim>
Tensor<2, dim>
FE_DG_DivFree<dim, spacedim>::shape_grad_grad(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const unsigned int i, const Point<dim> &p) const
{
  // Assert(this->is_primitive(), ExcShapeFunctionNotPrimitive());
  Assert(false, ExcNotImplemented());
  return Tensor<2, dim>();
}


template <int dim, int spacedim>
Tensor<2, dim>
FE_DG_DivFree<dim, spacedim>::shape_grad_grad_component(const unsigned int i, const Point<dim> &p, const unsigned int component) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<2, dim>();
}


template <int dim, int spacedim>
Tensor<2, dim>
FE_DG_DivFree<dim, spacedim>::shape_grad_grad_component(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const unsigned int i, const Point<dim> &p, const unsigned int component) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<2, dim>();
}


template <int dim, int spacedim>
void
FE_DG_DivFree<dim, spacedim>::fill_fe_values(const typename Triangulation<dim, spacedim>::cell_iterator & cell,
  const CellSimilarity::Similarity,
  const Quadrature<dim> &,
  const Mapping<dim, spacedim> &,
  const typename Mapping<dim, spacedim>::InternalDataBase &,
  const dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase        &fe_internal,
  dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &output_data) const
{
  Assert(fe_internal.update_each & update_quadrature_points, ExcInternalError());

  const unsigned int n_q_points = mapping_data.quadrature_points.size();
  double h = cell->diameter();
  Point<dim> c = cell->center();

  for (unsigned int i = 0; i < n_q_points; ++i)
  {
    // TODO: Use extent_in_direction() & mapping instead of diameter()
    const Point<dim> p = (Point<dim>) (mapping_data.quadrature_points[i] - c) / h;

    // Array size is dofs_per_cell * n_components
    for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
    {
      for (unsigned int j = 0; j < dim; ++j)
      {
        if (output_data.shape_function_to_row_table[k * dim + j] != numbers::invalid_unsigned_int)
        {
          if (fe_internal.update_each & update_values)
            output_data.shape_values[output_data.shape_function_to_row_table[k * dim + j]][i] = this->shape_value_component(cell, k, p, j, h);
          if (fe_internal.update_each & update_gradients)
            output_data.shape_gradients[output_data.shape_function_to_row_table[k * dim + j]][i] = this->shape_grad_component(cell, k, p, j);
        }
      }
    }
  }
}


template <int dim, int spacedim>
void
FE_DG_DivFree<dim, spacedim>::fill_fe_face_values(const typename Triangulation<dim, spacedim>::cell_iterator & cell,
  const unsigned int face_no,
  const Quadrature<dim - 1> &,
  const Mapping<dim, spacedim> &,
  const typename Mapping<dim, spacedim>::InternalDataBase &,
  const dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase        &fe_internal,
  dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &output_data) const
{
  Assert(fe_internal.update_each & update_quadrature_points, ExcInternalError());

  const unsigned int n_q_points = mapping_data.quadrature_points.size();
  double h = cell->diameter();
  Point<dim> c = cell->center();

  for (unsigned int i = 0; i < n_q_points; ++i)
  {
    // TODO: Use extent_in_direction() & mapping instead of diameter()
    const Point<dim> p = (Point<dim>) (mapping_data.quadrature_points[i] - c) / h;

    // Array size is dofs_per_cell * n_components
    for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
    {
      for (unsigned int j = 0; j < dim; ++j)
      {
        if (output_data.shape_function_to_row_table[k * dim + j] != numbers::invalid_unsigned_int)
        {
          if (fe_internal.update_each & update_values)
            output_data.shape_values[output_data.shape_function_to_row_table[k * dim + j]][i] = this->shape_value_component(cell, k, p, j, h);
          if (fe_internal.update_each & update_gradients)
            output_data.shape_gradients[output_data.shape_function_to_row_table[k * dim + j]][i] = this->shape_grad_component(cell, k, p, j);
        }
      }
    }
  }
}

template <int dim, int spacedim>
void FE_DG_DivFree<dim, spacedim>::fill_fe_subface_values(const typename Triangulation<dim, spacedim>::cell_iterator & cell,
  const unsigned int,
  const unsigned int,
  const Quadrature<dim - 1>                                             &,
  const Mapping<dim, spacedim> &,
  const typename Mapping<dim, spacedim>::InternalDataBase &,
  const dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase        &fe_internal,
  dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &output_data) const
{
  Assert(fe_internal.update_each & update_quadrature_points, ExcInternalError());

  const unsigned int n_q_points = mapping_data.quadrature_points.size();
  double h = cell->diameter();
  Point<dim> c = cell->center();

  for (unsigned int i = 0; i < n_q_points; ++i)
  {
    // TODO: Use extent_in_direction() & mapping instead of diameter()
    const Point<dim> p = (Point<dim>) (mapping_data.quadrature_points[i] - c) / h;

    // Array size is dofs_per_cell * n_components
    for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
    {
      for (unsigned int j = 0; j < dim; ++j)
      {
        if (output_data.shape_function_to_row_table[k * dim + j] != numbers::invalid_unsigned_int)
        {
          if (fe_internal.update_each & update_values)
            output_data.shape_values[output_data.shape_function_to_row_table[k * dim + j]][i] = this->shape_value_component(cell, k, p, j, h);
          if (fe_internal.update_each & update_gradients)
            output_data.shape_gradients[output_data.shape_function_to_row_table[k * dim + j]][i] = this->shape_grad_component(cell, k, p, j);
        }
      }
    }
  }
}


template <int dim, int spacedim>
void FE_DG_DivFree<dim, spacedim>::get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &x_source_fe, FullMatrix<double> &interpolation_matrix) const
{
  typedef              FiniteElement<dim, spacedim> FEE;
  AssertThrow((x_source_fe.get_name().find("FE_DG_DivFree<") == 0) || (dynamic_cast<const FE_DG_DivFree<dim, spacedim>*>(&x_source_fe) != 0), typename FEE::ExcInterpolationNotImplemented());

  Assert(interpolation_matrix.m() == 0, ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0, ExcDimensionMismatch(interpolation_matrix.n(), 0));
}


template <int dim, int spacedim>
void FE_DG_DivFree<dim, spacedim>::get_subface_interpolation_matrix(const FiniteElement<dim, spacedim> &x_source_fe,
  const unsigned int,
  FullMatrix<double>           &interpolation_matrix) const
{
  typedef              FiniteElement<dim, spacedim> FEE;
  AssertThrow((x_source_fe.get_name().find("FE_DG_DivFree<") == 0) || (dynamic_cast<const FE_DG_DivFree<dim, spacedim>*>(&x_source_fe) != 0), typename FEE::ExcInterpolationNotImplemented());

  Assert(interpolation_matrix.m() == 0, ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0, ExcDimensionMismatch(interpolation_matrix.n(), 0));
}


template <int dim, int spacedim>
bool FE_DG_DivFree<dim, spacedim>::hp_constraints_are_implemented() const
{
  return true;
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DG_DivFree<dim, spacedim>::hp_vertex_dof_identities(const FiniteElement<dim, spacedim> &fe_other) const
{
  // there are no such constraints for DGPNonparametric elements at all
  if (dynamic_cast<const FE_DG_DivFree<dim, spacedim>*>(&fe_other) != 0)
    return std::vector<std::pair<unsigned int, unsigned int> >();
  else
  {
    Assert(false, ExcNotImplemented());
    return std::vector<std::pair<unsigned int, unsigned int> >();
  }
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DG_DivFree<dim, spacedim>::hp_line_dof_identities(const FiniteElement<dim, spacedim> &fe_other) const
{
  // there are no such constraints for DGPNonparametric elements at all
  if (dynamic_cast<const FE_DG_DivFree<dim, spacedim>*>(&fe_other) != 0)
    return std::vector<std::pair<unsigned int, unsigned int> >();
  else
  {
    Assert(false, ExcNotImplemented());
    return std::vector<std::pair<unsigned int, unsigned int> >();
  }
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DG_DivFree<dim, spacedim>::hp_quad_dof_identities(const FiniteElement<dim, spacedim>        &fe_other) const
{
  // there are no such constraints for DGPNonparametric elements at all
  if (dynamic_cast<const FE_DG_DivFree<dim, spacedim>*>(&fe_other) != 0)
    return std::vector<std::pair<unsigned int, unsigned int> >();
  else
  {
    Assert(false, ExcNotImplemented());
    return std::vector<std::pair<unsigned int, unsigned int> >();
  }
}


template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_DG_DivFree<dim, spacedim>::compare_for_face_domination(const FiniteElement<dim, spacedim> &fe_other) const
{
  // check whether both are discontinuous elements, see the description of FiniteElementDomination::Domination
  if (dynamic_cast<const FE_DG_DivFree<dim, spacedim>*>(&fe_other) != 0)
    return FiniteElementDomination::no_requirements;

  Assert(false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}


template <int dim, int spacedim>
bool FE_DG_DivFree<dim, spacedim>::has_support_on_face(const unsigned int, const unsigned int) const
{
  return true;
}


template <int dim, int spacedim>
std::size_t FE_DG_DivFree<dim, spacedim>::memory_consumption() const
{
  Assert(false, ExcNotImplemented());
  return 0;
}

template <int dim, int spacedim>
bool FE_DG_DivFree<dim, spacedim>::is_constant(const unsigned int i) const
{
  return i < 3;
}


template <int dim, int spacedim>
unsigned int FE_DG_DivFree<dim, spacedim>::get_degree() const
{
  return this->degree;
}


template <int dim, int spacedim>
std::string FE_DG_DivFree<dim, spacedim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_DG_DivFree<" << Utilities::dim_string(dim, spacedim) << ">(" << this->degree << ")";
  return namebuf.str();
}

template class FE_DG_DivFree<3>;