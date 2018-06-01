// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include "feTaylor.h"

template <int dim, int spacedim>
FE_DG_Taylor<dim, spacedim>::FE_DG_Taylor(const unsigned int degree)
  :
  FiniteElement<dim, spacedim>(
    FiniteElementData<dim>(get_dpo_vector(degree), 1, degree,
      FiniteElementData<dim>::L2),
    std::vector<bool>(
      FiniteElementData<dim>(get_dpo_vector(degree), 1, degree).dofs_per_cell, true),
    std::vector<ComponentMask>(
      FiniteElementData<dim>(get_dpo_vector(degree), 1, degree).dofs_per_cell,
      std::vector<bool>(1, true))), FiniteElementIsConstantInterface<dim>(),
  polynomial_space(Polynomials::Monomial<double>::generate_complete_basis(degree))
{
  this->reinit_restriction_and_prolongation_matrices();
  // Fill prolongation matrices with embedding operators
  if (dim == spacedim)
  {
    FETools::compute_embedding_matrices(*this, this->prolongation);
    // Fill restriction matrices with L2-projection
    FETools::compute_projection_matrices(*this, this->restriction);
  }
}

template <int dim, int spacedim>
std::string FE_DG_Taylor<dim, spacedim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_DG_Taylor<" << Utilities::dim_string(dim, spacedim) << ">(" << this->degree << ")";
  return namebuf.str();
}


template <int dim, int spacedim>
FiniteElement<dim, spacedim> * FE_DG_Taylor<dim, spacedim>::clone() const
{
  return new FE_DG_Taylor<dim, spacedim>(*this);
}


template <int dim, int spacedim>
double FE_DG_Taylor<dim, spacedim>::shape_value(const unsigned int i, const Point<dim> &p) const
{
  Assert(false, ExcNotImplemented());
  return 0.;
}


template <int dim, int spacedim>
double
FE_DG_Taylor<dim, spacedim>::shape_value(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const unsigned int i, const Point<dim> &p) const
{
  Assert(false, ExcNotImplemented());
  return 0.;
}


template <int dim, int spacedim>
void
FE_DG_Taylor<dim, spacedim>::shape_values(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const std::vector< Point<dim> > &p, std::vector< std::vector<double> >& values) const
{
  Assert(false, ExcNotImplemented());
}


template <int dim, int spacedim>
double FE_DG_Taylor<dim, spacedim>::shape_value_component(const unsigned int i, const Point<dim> &p, const unsigned int component) const
{
  Assert(false, ExcNotImplemented());
  return 0.;
}


template <int dim, int spacedim>
double
FE_DG_Taylor<dim, spacedim>::shape_value_component(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const unsigned int i, const Point<dim> &p, const unsigned int component) const
{
  Assert(false, ExcNotImplemented());
  return 0.;
}


template <int dim, int spacedim>
Tensor<1, dim> FE_DG_Taylor<dim, spacedim>::shape_grad(const unsigned int i, const Point<dim> &p) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<1, dim>();
}


template <int dim, int spacedim>
Tensor<1, dim>
FE_DG_Taylor<dim, spacedim>::shape_grad(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const unsigned int i, const Point<dim> &p) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<1, dim>();
}


template <int dim, int spacedim>
Tensor<1, dim>
FE_DG_Taylor<dim, spacedim>::shape_grad_component(const unsigned int i, const Point<dim> &p, const unsigned int component) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<1, dim>();
}


template <int dim, int spacedim>
Tensor<1, dim>
FE_DG_Taylor<dim, spacedim>::shape_grad_component(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const unsigned int i, const Point<dim> &p, const unsigned int component) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<1, dim>();
}


template <int dim, int spacedim>
Tensor<2, dim>
FE_DG_Taylor<dim, spacedim>::shape_grad_grad(const unsigned int i, const Point<dim> &p) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<2, dim>();
}


template <int dim, int spacedim>
Tensor<2, dim>
FE_DG_Taylor<dim, spacedim>::shape_grad_grad(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const unsigned int i, const Point<dim> &p) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<2, dim>();
}


template <int dim, int spacedim>
Tensor<2, dim>
FE_DG_Taylor<dim, spacedim>::shape_grad_grad_component(const unsigned int i, const Point<dim> &p, const unsigned int component) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<2, dim>();
}


template <int dim, int spacedim>
Tensor<2, dim>
FE_DG_Taylor<dim, spacedim>::shape_grad_grad_component(const typename Triangulation<dim, spacedim>::cell_iterator & cell, const unsigned int i, const Point<dim> &p, const unsigned int component) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<2, dim>();
}


template <int dim, int spacedim>
std::vector<unsigned int>
FE_DG_Taylor<dim, spacedim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim + 1, static_cast<unsigned int>(0));
  dpo[dim] = deg + 1;
  for (unsigned int i = 1; i < dim; ++i)
  {
    dpo[dim] *= deg + 1 + i;
    dpo[dim] /= i + 1;
  }
  return dpo;
}


template <int dim, int spacedim>
UpdateFlags FE_DG_Taylor<dim, spacedim>::requires_update_flags(const UpdateFlags flags) const
{
  UpdateFlags out = flags;

  if (flags & (update_values | update_gradients | update_hessians))
    out |= update_quadrature_points;

  return out;
}


template <int dim, int spacedim>
typename FiniteElement<dim, spacedim>::InternalDataBase *
FE_DG_Taylor<dim, spacedim>::get_data(const UpdateFlags update_flags, const Mapping<dim, spacedim> &, const Quadrature<dim> &, dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &/*output_data*/) const
{
  typename FiniteElement<dim, spacedim>::InternalDataBase *data = new typename FiniteElement<dim, spacedim>::InternalDataBase;
  data->update_each = requires_update_flags(update_flags);
  return data;
}


template <int dim, int spacedim>
void
FE_DG_Taylor<dim, spacedim>::fill_fe_values(const typename Triangulation<dim, spacedim>::cell_iterator & cell,
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

  std::vector<double> values(fe_internal.update_each & update_values ? this->dofs_per_cell : 0);
  std::vector<Tensor<1, dim> > grads(fe_internal.update_each & update_gradients ? this->dofs_per_cell : 0);
  std::vector<Tensor<2, dim> > grad_grads(fe_internal.update_each & update_hessians ? this->dofs_per_cell : 0);
  std::vector<Tensor<3, dim> > empty_vector_of_3rd_order_tensors;//not used here, as well as elsewhere. not added everywhere!
  std::vector<Tensor<4, dim> > empty_vector_of_4th_order_tensors; //same as above and not well implemented in deal

  double h = cell->diameter();
  Point<dim> c = cell->center();

  if (fe_internal.update_each & (update_values | update_gradients))
    for (unsigned int i = 0; i < n_q_points; ++i)
    {
      const Point<dim> p = (Point<dim>)(mapping_data.quadrature_points[i] - c) / h;
      polynomial_space.compute(p, //mapping_data.quadrature_points[i],
        values, grads, grad_grads,
        empty_vector_of_3rd_order_tensors,
        empty_vector_of_4th_order_tensors);
      if (fe_internal.update_each & update_values)
        for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
          output_data.shape_values[k][i] = values[k];

      if (fe_internal.update_each & update_gradients)
        for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
          output_data.shape_gradients[k][i] = grads[k] / h;

      if (fe_internal.update_each & update_hessians)
        for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
          output_data.shape_hessians[k][i] = grad_grads[k] / h / h;
    }
}


template <int dim, int spacedim>
void
FE_DG_Taylor<dim, spacedim>::fill_fe_face_values(const typename Triangulation<dim, spacedim>::cell_iterator & cell,
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

  std::vector<double> values(fe_internal.update_each & update_values ? this->dofs_per_cell : 0);
  std::vector<Tensor<1, dim> > grads(fe_internal.update_each & update_gradients ? this->dofs_per_cell : 0);
  std::vector<Tensor<2, dim> > grad_grads(fe_internal.update_each & update_hessians ? this->dofs_per_cell : 0);
  std::vector<Tensor<3, dim> > empty_vector_of_3rd_order_tensors;
  std::vector<Tensor<4, dim> > empty_vector_of_4th_order_tensors;

  double h = cell->diameter();
  Point<dim> c = cell->center();

  if (fe_internal.update_each & (update_values | update_gradients))
  {
    for (unsigned int i = 0; i < n_q_points; ++i)
    {
      const Point<dim> p = (Point<dim>) (mapping_data.quadrature_points[i] - c) / h;
      polynomial_space.compute(p, //mapping_data.quadrature_points[i],
        values, grads, grad_grads,
        empty_vector_of_3rd_order_tensors,
        empty_vector_of_4th_order_tensors);
      if (fe_internal.update_each & update_values)
        for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
          output_data.shape_values[k][i] = values[k];

      if (fe_internal.update_each & update_gradients)
        for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
          output_data.shape_gradients[k][i] = grads[k] / h;

      if (fe_internal.update_each & update_hessians)
        for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
          output_data.shape_hessians[k][i] = grad_grads[k] / h / h;
    }
  }
}


template <int dim, int spacedim>
void FE_DG_Taylor<dim, spacedim>::fill_fe_subface_values(const typename Triangulation<dim, spacedim>::cell_iterator & cell,
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

  std::vector<double> values(fe_internal.update_each & update_values ? this->dofs_per_cell : 0);
  std::vector<Tensor<1, dim> > grads(fe_internal.update_each & update_gradients ? this->dofs_per_cell : 0);
  std::vector<Tensor<2, dim> > grad_grads(fe_internal.update_each & update_hessians ? this->dofs_per_cell : 0);
  std::vector<Tensor<3, dim> > empty_vector_of_3rd_order_tensors;
  std::vector<Tensor<4, dim> > empty_vector_of_4th_order_tensors;

  double h = cell->diameter();

  if (fe_internal.update_each & (update_values | update_gradients))
    for (unsigned int i = 0; i < n_q_points; ++i)
    {
      const Point<dim> p = (Point<dim>)(mapping_data.quadrature_points[i] - cell->center()) / h;
      polynomial_space.compute(p, //mapping_data.quadrature_points[i],
        values, grads, grad_grads,
        empty_vector_of_3rd_order_tensors,
        empty_vector_of_4th_order_tensors);
      if (fe_internal.update_each & update_values)
        for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
          output_data.shape_values[k][i] = values[k];

      if (fe_internal.update_each & update_gradients)
        for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
          output_data.shape_gradients[k][i] = grads[k] / h;

      if (fe_internal.update_each & update_hessians)
        for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
          output_data.shape_hessians[k][i] = grad_grads[k] / h / h;
    }
}


template <int dim, int spacedim>
void FE_DG_Taylor<dim, spacedim>::get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &x_source_fe, FullMatrix<double> &interpolation_matrix) const
{
  typedef              FiniteElement<dim, spacedim> FEE;
  AssertThrow((x_source_fe.get_name().find("FE_DG_Taylor<") == 0) || (dynamic_cast<const FE_DG_Taylor<dim, spacedim>*>(&x_source_fe) != 0), typename FEE::ExcInterpolationNotImplemented());

  Assert(interpolation_matrix.m() == 0, ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0, ExcDimensionMismatch(interpolation_matrix.n(), 0));
}


template <int dim, int spacedim>
void FE_DG_Taylor<dim, spacedim>::get_subface_interpolation_matrix(const FiniteElement<dim, spacedim> &x_source_fe,
  const unsigned int,
  FullMatrix<double>           &interpolation_matrix) const
{
  typedef              FiniteElement<dim, spacedim> FEE;
  AssertThrow((x_source_fe.get_name().find("FE_DG_Taylor<") == 0) || (dynamic_cast<const FE_DG_Taylor<dim, spacedim>*>(&x_source_fe) != 0), typename FEE::ExcInterpolationNotImplemented());

  Assert(interpolation_matrix.m() == 0, ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0, ExcDimensionMismatch(interpolation_matrix.n(), 0));
}


template <int dim, int spacedim>
bool FE_DG_Taylor<dim, spacedim>::hp_constraints_are_implemented() const
{
  return true;
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DG_Taylor<dim, spacedim>::hp_vertex_dof_identities(const FiniteElement<dim, spacedim> &fe_other) const
{
  // there are no such constraints for DGPNonparametric elements at all
  if (dynamic_cast<const FE_DG_Taylor<dim, spacedim>*>(&fe_other) != 0)
    return std::vector<std::pair<unsigned int, unsigned int> >();
  else
  {
    Assert(false, ExcNotImplemented());
    return std::vector<std::pair<unsigned int, unsigned int> >();
  }
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DG_Taylor<dim, spacedim>::hp_line_dof_identities(const FiniteElement<dim, spacedim> &fe_other) const
{
  // there are no such constraints for DGPNonparametric elements at all
  if (dynamic_cast<const FE_DG_Taylor<dim, spacedim>*>(&fe_other) != 0)
    return std::vector<std::pair<unsigned int, unsigned int> >();
  else
  {
    Assert(false, ExcNotImplemented());
    return std::vector<std::pair<unsigned int, unsigned int> >();
  }
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DG_Taylor<dim, spacedim>::hp_quad_dof_identities(const FiniteElement<dim, spacedim>        &fe_other) const
{
  // there are no such constraints for DGPNonparametric elements at all
  if (dynamic_cast<const FE_DG_Taylor<dim, spacedim>*>(&fe_other) != 0)
    return std::vector<std::pair<unsigned int, unsigned int> >();
  else
  {
    Assert(false, ExcNotImplemented());
    return std::vector<std::pair<unsigned int, unsigned int> >();
  }
}


template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_DG_Taylor<dim, spacedim>::compare_for_face_domination(const FiniteElement<dim, spacedim> &fe_other) const
{
  // check whether both are discontinuous elements, see the description of FiniteElementDomination::Domination
  if (dynamic_cast<const FE_DG_Taylor<dim, spacedim>*>(&fe_other) != 0)
    return FiniteElementDomination::no_requirements;

  Assert(false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}


template <int dim, int spacedim>
bool FE_DG_Taylor<dim, spacedim>::has_support_on_face(const unsigned int shape_index, const unsigned int face_index) const
{
  return true;
}

template <int dim, int spacedim>
bool FE_DG_Taylor<dim, spacedim>::is_constant(const unsigned int i) const
{
  return i < 1;
}


template <int dim, int spacedim>
std::size_t FE_DG_Taylor<dim, spacedim>::memory_consumption() const
{
  Assert(false, ExcNotImplemented());
  return 0;
}


template <int dim, int spacedim>
unsigned int FE_DG_Taylor<dim, spacedim>::get_degree() const
{
  return this->degree;
}

template class FE_DG_Taylor<3>;