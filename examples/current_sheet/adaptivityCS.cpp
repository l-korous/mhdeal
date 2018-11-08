#include "adaptivityCS.h"

template <int dim>
AdaptivityCS<dim>::AdaptivityCS(Parameters<dim>& parameters, MPI_Comm& mpi_communicator) :
  Adaptivity<dim>(parameters, mpi_communicator),
  last_time_step(0),
  adaptivity_step(0),
  mag(dim + 2)
{
}

template <int dim>
void AdaptivityCS<dim>::calculate_jumps(TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler, const Mapping<dim>& mapping, Vector<double>& gradient_indicator)
{
  FEValuesExtractors::Scalar scalars[dim];
  scalars[0].component = 0;
  scalars[1].component = 4;

  const QGauss<dim - 1> face_quadrature(1);
  UpdateFlags face_update_flags = UpdateFlags(update_values | update_JxW_values | update_gradients);
  FEFaceValues<dim> fe_v_face(mapping, dof_handler.get_fe(), face_quadrature, face_update_flags);
  FESubfaceValues<dim> fe_v_subface(mapping, dof_handler.get_fe(), face_quadrature, face_update_flags);
  FEFaceValues<dim> fe_v_face_neighbor(mapping, dof_handler.get_fe(), face_quadrature, update_values | update_gradients);
  int n_quadrature_points_face = face_quadrature.get_points().size();
  int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
  this->dof_indices.resize(dofs_per_cell);
  this->dof_indices_neighbor.resize(dofs_per_cell);
  for (unsigned int i = 0; i < dofs_per_cell; ++i)
  {
    const unsigned int component_i = dof_handler.get_fe().system_to_base_index(i).first.first;
    is_primitive[i] = dof_handler.get_fe().is_primitive(i);
    if (is_primitive[i])
      component_ii[i] = dof_handler.get_fe().system_to_component_index(i).first;
    else
      component_ii[i] = 999;
  }

  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
  {
    if (!cell->is_locally_owned())
      continue;
    Point<dim> jump;
    Point<dim> area;
    cell->get_dof_indices(this->dof_indices);
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
    {
      typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
      if (!face->at_boundary())
      {
        Assert(cell->neighbor(face_no).state() == IteratorState::valid, ExcInternalError());
        typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
        std::vector<double> u(n_quadrature_points_face);
        std::vector<double> u_neighbor(n_quadrature_points_face);
        std::vector<std::array<std::array<double, dim>, dim> > u_mag;
        std::vector<std::array<std::array<double, dim>, dim> > u_neighbor_mag;
        u_mag.resize(n_quadrature_points_face);
        u_neighbor_mag.resize(n_quadrature_points_face);
        if (face->has_children())
        {
          unsigned int neighbor2 = cell->neighbor_face_no(face_no);
          for (unsigned int subface_no = 0; subface_no < face->number_of_children(); ++subface_no)
          {
            typename DoFHandler<dim>::cell_iterator neighbor_child = cell->neighbor_child_on_subface(face_no, subface_no);
            Assert(!neighbor_child->has_children(), ExcInternalError());
            fe_v_subface.reinit(cell, face_no, subface_no);
            fe_v_face_neighbor.reinit(neighbor_child, neighbor2);
            neighbor_child->get_dof_indices(dof_indices_neighbor);
            const std::vector<double> &JxW = fe_v_subface.get_JxW_values();
            for (unsigned int x = 0; x < n_quadrature_points_face; ++x)
              area[face_no / 2] += JxW[x];
            if (this->parameters.polynomial_order_dg == 0)
            {
              for (int scalar_i = 0; scalar_i < 2; scalar_i++)
              {
                fe_v_subface[scalars[scalar_i]].get_function_values(solution, u);
                fe_v_face_neighbor[scalars[scalar_i]].get_function_values(solution, u_neighbor);

                for (unsigned int x = 0; x < n_quadrature_points_face; ++x)
                  jump[face_no / 2] += std::fabs(u[x] - u_neighbor[x]) * JxW[x];
              }
            }
            else
            {
              for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
                for (int d = 0; d < dim; d++)
                  for (int e = 0; e < dim; e++)
                    u_mag[q][d][e] = u_neighbor_mag[q][d][e] = 0.;

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
                {
                  // Plus
                  if (!is_primitive[i])
                  {
                    Tensor<2, dim> fe_v_grad = fe_v_subface[mag].gradient(i, q);
                    for (int d = 0; d < dim; d++)
                      for (int e = 0; e < dim; e++)
                        u_mag[q][d][e] += solution(this->dof_indices[i]) * fe_v_grad[d][e];
                  }
                  else if(component_ii[i] >= 5)
                    for (int d = 0; d < dim; d++)
                      u_mag[q][component_ii[i] - 5][d] += solution(this->dof_indices[i]) * fe_v_subface.shape_grad(i, q)[d];

                  // Minus
                  if (!is_primitive[i])
                  {
                    Tensor<2, dim> fe_v_grad = fe_v_face_neighbor[mag].gradient(i, q);
                    for (int d = 0; d < dim; d++)
                      for (int e = 0; e < dim; e++)
                        u_neighbor_mag[q][d][e] += solution(this->dof_indices_neighbor[i]) * fe_v_grad[d][e];
                  }
                  else if (component_ii[i] >= 5)
                    for (int d = 0; d < dim; d++)
                      u_neighbor_mag[q][component_ii[i] - 5][d] += solution(this->dof_indices_neighbor[i]) * fe_v_face_neighbor.shape_grad(i, q)[d];
                }
              }

              for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
              {
                std::array<double, dim> curl_ = { u_mag[q][2][1] - u_mag[q][1][2], u_mag[q][0][2] - u_mag[q][2][0], u_mag[q][1][0] - u_mag[q][0][1] };
                u[q] = curl_[0] * curl_[0] + curl_[1] * curl_[1] + curl_[2] * curl_[2];

                std::array<double, dim> curl_neighbor_ = { u_neighbor_mag[q][2][1] - u_neighbor_mag[q][1][2], u_neighbor_mag[q][0][2] - u_neighbor_mag[q][2][0], u_neighbor_mag[q][1][0] - u_neighbor_mag[q][0][1] };
                u_neighbor[q] = curl_neighbor_[0] * curl_neighbor_[0] + curl_neighbor_[1] * curl_neighbor_[1] + curl_neighbor_[2] * curl_neighbor_[2];

                jump[face_no / 2] += std::fabs(u[q] - u_neighbor[q]) * JxW[q];
              }
            }
          }
        }
        else
        {
          if (!cell->neighbor_is_coarser(face_no))
          {
            unsigned int neighbor2 = cell->neighbor_of_neighbor(face_no);
            fe_v_face.reinit(cell, face_no);
            fe_v_face_neighbor.reinit(neighbor, neighbor2);
            neighbor->get_dof_indices(dof_indices_neighbor);
            const std::vector<double> &JxW = fe_v_face.get_JxW_values();
            for (unsigned int x = 0; x < n_quadrature_points_face; ++x)
              area[face_no / 2] += JxW[x];
            if (this->parameters.polynomial_order_dg == 0)
            {
              for (int scalar_i = 0; scalar_i < 2; scalar_i++)
              {
                fe_v_face[scalars[scalar_i]].get_function_values(solution, u);
                fe_v_face_neighbor[scalars[scalar_i]].get_function_values(solution, u_neighbor);

                for (unsigned int x = 0; x < n_quadrature_points_face; ++x)
                  jump[face_no / 2] += std::fabs(u[x] - u_neighbor[x]) * JxW[x];
              }
            }
            else
            {
              for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
                for (int d = 0; d < dim; d++)
                  for (int e = 0; e < dim; e++)
                    u_mag[q][d][e] = u_neighbor_mag[q][d][e] = 0.;

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
                {
                  // Plus
                  if (!is_primitive[i])
                  {
                    Tensor<2, dim> fe_v_grad = fe_v_face[mag].gradient(i, q);
                    for (int d = 0; d < dim; d++)
                      for (int e = 0; e < dim; e++)
                        u_mag[q][d][e] += solution(this->dof_indices[i]) * fe_v_grad[d][e];
                  }
                  else if (component_ii[i] >= 5)
                    for (int d = 0; d < dim; d++)
                      u_mag[q][component_ii[i] - 5][d] += solution(this->dof_indices[i]) * fe_v_face.shape_grad(i, q)[d];

                  // Minus
                  if (!is_primitive[i])
                  {
                    Tensor<2, dim> fe_v_grad = fe_v_face_neighbor[mag].gradient(i, q);
                    for (int d = 0; d < dim; d++)
                      for (int e = 0; e < dim; e++)
                        u_neighbor_mag[q][d][e] += solution(this->dof_indices_neighbor[i]) * fe_v_grad[d][e];
                  }
                  else if (component_ii[i] >= 5)
                    for (int d = 0; d < dim; d++)
                      u_neighbor_mag[q][component_ii[i] - 5][d] += solution(this->dof_indices_neighbor[i]) * fe_v_face_neighbor.shape_grad(i, q)[d];
                }
              }

              for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
              {
                std::array<double, dim> curl_ = { u_mag[q][2][1] - u_mag[q][1][2], u_mag[q][0][2] - u_mag[q][2][0], u_mag[q][1][0] - u_mag[q][0][1] };
                u[q] = curl_[0] * curl_[0] + curl_[1] * curl_[1] + curl_[2] * curl_[2];

                std::array<double, dim> curl_neighbor_ = { u_neighbor_mag[q][2][1] - u_neighbor_mag[q][1][2], u_neighbor_mag[q][0][2] - u_neighbor_mag[q][2][0], u_neighbor_mag[q][1][0] - u_neighbor_mag[q][0][1] };
                u_neighbor[q] = curl_neighbor_[0] * curl_neighbor_[0] + curl_neighbor_[1] * curl_neighbor_[1] + curl_neighbor_[2] * curl_neighbor_[2];

                jump[face_no / 2] += std::fabs(u[q] - u_neighbor[q]) * JxW[q];
              }
            }
          }
          else //i.e. neighbor is coarser than cell
          {
            std::pair<unsigned int, unsigned int> neighbor_face_subface
              = cell->neighbor_of_coarser_neighbor(face_no);
            Assert(neighbor_face_subface.first < GeometryInfo<dim>::faces_per_cell, ExcInternalError());
            Assert(neighbor_face_subface.second < neighbor->face(neighbor_face_subface.first)->number_of_children(),
              ExcInternalError());
            Assert(neighbor->neighbor_child_on_subface(neighbor_face_subface.first, neighbor_face_subface.second)
              == cell, ExcInternalError());
            fe_v_face.reinit(cell, face_no);
            fe_v_subface.reinit(neighbor, neighbor_face_subface.first,
              neighbor_face_subface.second);
            neighbor->get_dof_indices(dof_indices_neighbor);
            const std::vector<double> &JxW = fe_v_face.get_JxW_values();
            for (unsigned int x = 0; x < n_quadrature_points_face; ++x)
              area[face_no / 2] += JxW[x];
            if (this->parameters.polynomial_order_dg == 0)
            {
            for (int scalar_i = 0; scalar_i < 2; scalar_i++)
            {
              fe_v_face[scalars[scalar_i]].get_function_values(solution, u);
              fe_v_subface[scalars[scalar_i]].get_function_values(solution, u_neighbor);

              for (unsigned int x = 0; x < n_quadrature_points_face; ++x)
                jump[face_no / 2] += std::fabs(u[x] - u_neighbor[x]) * JxW[x];
            }
            }
            else
            {
              for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
                for (int d = 0; d < dim; d++)
                  for (int e = 0; e < dim; e++)
                    u_mag[q][d][e] = u_neighbor_mag[q][d][e] = 0.;

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
                {
                  // Plus
                  if (!is_primitive[i])
                  {
                    Tensor<2, dim> fe_v_grad = fe_v_face[mag].gradient(i, q);
                    for (int d = 0; d < dim; d++)
                      for (int e = 0; e < dim; e++)
                        u_mag[q][d][e] += solution(this->dof_indices[i]) * fe_v_grad[d][e];
                  }
                  else if (component_ii[i] >= 5)
                    for (int d = 0; d < dim; d++)
                      u_mag[q][component_ii[i] - 5][d] += solution(this->dof_indices[i]) * fe_v_face.shape_grad(i, q)[d];

                  // Minus
                  if (!is_primitive[i])
                  {
                    Tensor<2, dim> fe_v_grad = fe_v_subface[mag].gradient(i, q);
                    for (int d = 0; d < dim; d++)
                      for (int e = 0; e < dim; e++)
                        u_neighbor_mag[q][d][e] += solution(this->dof_indices_neighbor[i]) * fe_v_grad[d][e];
                  }
                  else if (component_ii[i] >= 5)
                    for (int d = 0; d < dim; d++)
                      u_neighbor_mag[q][component_ii[i] - 5][d] += solution(this->dof_indices_neighbor[i]) * fe_v_subface.shape_grad(i, q)[d];
                }
              }

              for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
              {
                std::array<double, dim> curl_ = { u_mag[q][2][1] - u_mag[q][1][2], u_mag[q][0][2] - u_mag[q][2][0], u_mag[q][1][0] - u_mag[q][0][1] };
                u[q] = curl_[0] * curl_[0] + curl_[1] * curl_[1] + curl_[2] * curl_[2];

                std::array<double, dim> curl_neighbor_ = { u_neighbor_mag[q][2][1] - u_neighbor_mag[q][1][2], u_neighbor_mag[q][0][2] - u_neighbor_mag[q][2][0], u_neighbor_mag[q][1][0] - u_neighbor_mag[q][0][1] };
                u_neighbor[q] = curl_neighbor_[0] * curl_neighbor_[0] + curl_neighbor_[1] * curl_neighbor_[1] + curl_neighbor_[2] * curl_neighbor_[2];

                jump[face_no / 2] += std::fabs(u[q] - u_neighbor[q]) * JxW[q];
              }
            }
          }
        }
      }
    }
    double average_jumps[dim];
    double sum_of_average_jumps = 0.;
    for (unsigned int i = 0; i < dim; ++i)
    {
      average_jumps[i] = jump(i) / area(i);
      sum_of_average_jumps += average_jumps[i];
    }
    for (int i = 0; i < this->parameters.volume_factor; i++)
      sum_of_average_jumps *= cell->diameter();
    gradient_indicator(cell->active_cell_index()) = sum_of_average_jumps;
  }
}

template <int dim>
bool AdaptivityCS<dim>::refine_mesh(int time_step, double time, TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler,
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& triangulation
#else
  Triangulation<dim>& triangulation
#endif
  , const Mapping<dim>& mapping)
{
  if (time_step % this->parameters.refine_every_nth_time_step)
    return false;
  if (++adaptivity_step > (time_step == 0 ? this->parameters.perform_n_initial_refinements : 1))
  {
    adaptivity_step = 0;
    return false;
  }
  Vector<double> gradient_indicator(triangulation.n_active_cells());
  calculate_jumps(solution, dof_handler, mapping, gradient_indicator);

  int max_calls_ = this->parameters.max_cells + (int)std::floor(time * this->parameters.max_cells * this->parameters.time_interval_max_cells_multiplicator / this->parameters.final_time);
  GridRefinement::refine_and_coarsen_fixed_fraction(triangulation, gradient_indicator, this->parameters.refine_threshold, this->parameters.coarsen_threshold, max_calls_);

  triangulation.prepare_coarsening_and_refinement();

  return true;
}

template class AdaptivityCS<3>;
