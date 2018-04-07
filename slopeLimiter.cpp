#include "slopeLimiter.h"
#include "problem.h"

template <EquationsType equationsType, int dim>
void VertexBasedSlopeLimiter<equationsType, dim>::postprocess(TrilinosWrappers::MPI::Vector& current_limited_solution, TrilinosWrappers::MPI::Vector& current_unlimited_solution)
{
  int cell_count = 0;
  // Loop through all cells.
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
  {
    if (!cell->is_locally_owned())
      continue;

    bool u_c_set[Equations<equationsType, dim>::n_components];
    for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
      u_c_set[i] = false;

    double u_c[Equations<equationsType, dim>::n_components];
    cell->get_dof_indices(dof_indices);

    PostprocessData* data = 0;
    auto it = this->postprocessData.find(cell->active_cell_index());
    if (it != this->postprocessData.end())
      data = &(it->second);
    else
    {
      data = &(((postprocessData.insert(std::make_pair(cell->active_cell_index(), PostprocessData()))).first)->second);
      for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
        u_c_set[i] = false;
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        if (!is_primitive[i])
          data->lambda_indices_to_multiply_all_B_components.push_back(dof_indices[i]);
        else
        {
          // Here we rely on the fact, that the constant basis fn is the first one and all other basis fns come after.
          if (!u_c_set[component_ii[i]])
            u_c_set[component_ii[i]] = true;
          else
            data->lambda_indices_to_multiply[component_ii[i]].push_back(dof_indices[i]);
        }
      }

      data->center = cell->center();
      for (unsigned int vertex_i = 0; vertex_i < GeometryInfo<dim>::vertices_per_cell; ++vertex_i)
      {
        data->vertexIndex[vertex_i] = cell->vertex_index(vertex_i);
        data->vertexPoint[vertex_i] = data->center + (1. - NEGLIGIBLE) * (cell->vertex(vertex_i) - data->center);

        unsigned short neighbor_i = 0;

        for (auto neighbor_element : GridTools::find_cells_adjacent_to_vertex(triangulation, data->vertexIndex[vertex_i]))
        {
          typename DoFHandler<dim>::active_cell_iterator neighbor(&triangulation, neighbor_element->level(), neighbor_element->index(), &dof_handler);
          if (neighbor->active_cell_index() != cell->active_cell_index())
          {
            data->neighbor_dof_indices[vertex_i][neighbor_i].resize(dofs_per_cell);
            neighbor->get_dof_indices(data->neighbor_dof_indices[vertex_i][neighbor_i++]);
          }
        }

        // If on boundary which is periodic, get also the proper periodic neighbors
        for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
        {
          if (cell->at_boundary(face_no) && Problem<equationsType, dim>::is_periodic_boundary(cell->face(face_no)->boundary_id(), this->parameters))
          {
            TriaIterator<TriaAccessor<dim - 1, dim, dim> > face = cell->face(face_no);
            for (unsigned int face_i = 0; face_i < GeometryInfo<dim>::vertices_per_face; ++face_i)
            {
              // Only now we know this vertex is at a periodic boundary
              if (face->vertex_index(face_i) == data->vertexIndex[vertex_i])
              {
                const DealIIExtensions::FacePair<dim>&  face_pair = periodic_cell_map.find(std::make_pair(cell, face_no))->second;
                data->neighbor_dof_indices[vertex_i][neighbor_i].resize(dofs_per_cell);
                typename DoFHandler<dim>::active_cell_iterator neighbor(cell);
                auto this_cell_index = cell->active_cell_index();
                auto zeroth_found_cell_index = (*(face_pair.cell[0])).active_cell_index();
                neighbor = ((zeroth_found_cell_index == this_cell_index && face_no == face_pair.face_idx[0]) ? face_pair.cell[1] : face_pair.cell[0]);
                neighbor->get_dof_indices(data->neighbor_dof_indices[vertex_i][neighbor_i++]);
              }
            }
          }
        }
      }
    }

    // Cell center value we must find in any case (new data or reused)
    // - let us reuse this array for that.
    for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
      u_c_set[i] = false;
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      if (is_primitive[i])
      {
        // Here we rely on the fact, that the constant basis fn is the first one and that all other basis fns have zero mean.
        if (!u_c_set[component_ii[i]])
        {
          u_c[component_ii[i]] = current_unlimited_solution(dof_indices[i]);
          u_c_set[component_ii[i]] = true;
        }
      }
    }

    if (parameters.debug)
      std::cout << "cell: " << ++cell_count << " - center: " << data->center << ", values: " << u_c[0] << ", " << u_c[1] << ", " << u_c[2] << ", " << u_c[3] << ", " << u_c[4] << std::endl;

    double alpha_e[Equations<equationsType, dim>::n_components];
    for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
      alpha_e[i] = 1.;
    for (unsigned int vertex_i = 0; vertex_i < GeometryInfo<dim>::vertices_per_cell; ++vertex_i)
    {
      // (!!!) Find out u_i
      Vector<double> u_i(Equations<equationsType, dim>::n_components);
      const Point<dim> p_cell = mapping.transform_real_to_unit_cell(cell, data->vertexPoint[vertex_i]);
      const Quadrature<dim> one_point_quadrature(GeometryInfo<dim>::project_to_unit_cell(p_cell));
      FEValues<dim> fe_values(mapping, fe, one_point_quadrature, update_values);
      fe_values.reinit(cell);
      std::vector<Vector<double> > u_value(1, Vector<double>(fe.n_components()));
      fe_values.get_function_values(current_unlimited_solution, u_value);
      u_i = u_value[0];

      if (parameters.debug)
      {
        std::cout << "\tv_i: " << cell->vertex(vertex_i) << ", values: ";
        for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
          std::cout << u_i[i] << (i == Equations<equationsType, dim>::n_components - 1 ? "" : ", ");
        std::cout << std::endl;
      }

      // Init u_i_min, u_i_max
      double u_i_min[Equations<equationsType, dim>::n_components];
      double u_i_max[Equations<equationsType, dim>::n_components];
      for (int k = 0; k < Equations<equationsType, dim>::n_components; k++)
      {
        u_i_min[k] = u_c[k];
        u_i_max[k] = u_c[k];
      }

      // For all vertices -> v_i
      for (auto dof_indices_neighbor : data->neighbor_dof_indices[vertex_i])
      {
        if (dof_indices_neighbor.size() == 0)
          continue;
        bool u_i_extrema_set[Equations<equationsType, dim>::n_components];
        for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
          u_i_extrema_set[i] = false;

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          if (is_primitive[i])
          {
            // Here we rely on the fact, that the constant basis fn is the first one.
            if (!u_i_extrema_set[component_ii[i]])
            {
              double val = current_unlimited_solution(dof_indices_neighbor[i]);
              if (parameters.debug)
              {
                if (val < u_i_min[component_ii[i]])
                  std::cout << "\tdecreasing u_i_min to: " << val << std::endl;
                if (val > u_i_max[component_ii[i]])
                  std::cout << "\tincreasing u_i_max to: " << val << std::endl;
              }
              u_i_min[component_ii[i]] = std::min(u_i_min[component_ii[i]], val);
              u_i_max[component_ii[i]] = std::max(u_i_max[component_ii[i]], val);
              u_i_extrema_set[component_ii[i]] = true;
            }
          }
        }
      }

      // Based on u_i_min, u_i_max, u_i, get alpha_e
      for (int k = 0; k < Equations<equationsType, dim>::n_components; k++)
        if (std::abs((u_c[k] - u_i[k]) / u_c[k]) > NEGLIGIBLE)
        {
          alpha_e[k] = std::min(alpha_e[k], ((u_i[k] - u_c[k]) > 0.) ? std::min(1.0, (u_i_max[k] - u_c[k]) / (u_i[k] - u_c[k])) : std::min(1.0, (u_i_min[k] - u_c[k]) / (u_i[k] - u_c[k])));
          if (parameters.debug)
            std::cout << "\talpha_e[" << k << "]: " << alpha_e[k] << std::endl;
        }
    }

    for (int k = 0; k < Equations<equationsType, dim>::n_components; k++)
      for (int i = 0; i < data->lambda_indices_to_multiply[k].size(); i++)
        current_limited_solution(data->lambda_indices_to_multiply[k][i]) *= alpha_e[k];

    double alpha_e_B = std::min(std::min(alpha_e[5], alpha_e[6]), alpha_e[7]);
    for (int i = 0; i < data->lambda_indices_to_multiply_all_B_components.size(); i++)
      current_limited_solution(data->lambda_indices_to_multiply_all_B_components[i]) *= alpha_e_B;
  }
}

template <EquationsType equationsType, int dim>
void BarthJespersenSlopeLimiter<equationsType, dim>::postprocess(TrilinosWrappers::MPI::Vector& current_limited_solution, TrilinosWrappers::MPI::Vector& current_unlimited_solution)
{
  int cell_count = 0;
  // Loop through all cells.
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
  {
    if (!cell->is_locally_owned())
      continue;

    bool u_c_set[Equations<equationsType, dim>::n_components];
    for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
      u_c_set[i] = false;

    double u_c[Equations<equationsType, dim>::n_components];
    cell->get_dof_indices(dof_indices);

    PostprocessData* data = 0;
    auto it = this->postprocessData.find(cell->active_cell_index());
    if (it != this->postprocessData.end())
      data = &(it->second);
    else
    {
      data = &(((postprocessData.insert(std::make_pair(cell->active_cell_index(), PostprocessData()))).first)->second);
      for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
        u_c_set[i] = false;
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        if (!is_primitive[i])
          data->lambda_indices_to_multiply_all_B_components.push_back(dof_indices[i]);
        else
        {
          // Here we rely on the fact, that the constant basis fn is the first one and all other basis fns come after.
          if (!u_c_set[component_ii[i]])
            u_c_set[component_ii[i]] = true;
          else
            data->lambda_indices_to_multiply[component_ii[i]].push_back(dof_indices[i]);
        }
      }

      data->center = cell->center();
      for (unsigned int vertex_i = 0; vertex_i < GeometryInfo<dim>::vertices_per_cell; ++vertex_i)
      {
        data->vertexIndex[vertex_i] = cell->vertex_index(vertex_i);
        data->vertexPoint[vertex_i] = data->center + (1. - NEGLIGIBLE) * (cell->vertex(vertex_i) - data->center);

        unsigned short neighbor_i = 0;

        for (auto neighbor_element : GridTools::find_cells_adjacent_to_vertex(triangulation, data->vertexIndex[vertex_i]))
        {
          typename DoFHandler<dim>::active_cell_iterator neighbor(&triangulation, neighbor_element->level(), neighbor_element->index(), &dof_handler);
          if (neighbor->active_cell_index() != cell->active_cell_index())
          {
            data->neighbor_dof_indices[vertex_i][neighbor_i].resize(dofs_per_cell);
            neighbor->get_dof_indices(data->neighbor_dof_indices[vertex_i][neighbor_i++]);
          }
        }

        // If on boundary which is periodic, get also the proper periodic neighbors
        for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
        {
          if (cell->at_boundary(face_no) && Problem<equationsType, dim>::is_periodic_boundary(cell->face(face_no)->boundary_id(), this->parameters))
          {
            TriaIterator<TriaAccessor<dim - 1, dim, dim> > face = cell->face(face_no);
            for (unsigned int face_i = 0; face_i < GeometryInfo<dim>::vertices_per_face; ++face_i)
            {
              // Only now we know this vertex is at a periodic boundary
              if (face->vertex_index(face_i) == data->vertexIndex[vertex_i])
              {
                const DealIIExtensions::FacePair<dim>&  face_pair = periodic_cell_map.find(std::make_pair(cell, face_no))->second;
                data->neighbor_dof_indices[vertex_i][neighbor_i].resize(dofs_per_cell);
                typename DoFHandler<dim>::active_cell_iterator neighbor(cell);
                auto this_cell_index = cell->active_cell_index();
                auto zeroth_found_cell_index = (*(face_pair.cell[0])).active_cell_index();
                neighbor = ((zeroth_found_cell_index == this_cell_index && face_no == face_pair.face_idx[0]) ? face_pair.cell[1] : face_pair.cell[0]);
                neighbor->get_dof_indices(data->neighbor_dof_indices[vertex_i][neighbor_i++]);
              }
            }
          }
        }
      }
    }

    // Cell center value we must find in any case (new data or reused)
    // - let us reuse this array for that.
    for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
      u_c_set[i] = false;
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      if (is_primitive[i])
      {
        // Here we rely on the fact, that the constant basis fn is the first one and that all other basis fns have zero mean.
        if (!u_c_set[component_ii[i]])
        {
          u_c[component_ii[i]] = current_unlimited_solution(dof_indices[i]);
          u_c_set[component_ii[i]] = true;
        }
      }
    }

    if (parameters.debug)
      std::cout << "cell: " << ++cell_count << " - center: " << data->center << ", values: " << u_c[0] << ", " << u_c[1] << ", " << u_c[2] << ", " << u_c[3] << ", " << u_c[4] << std::endl;

    double alpha_e[Equations<equationsType, dim>::n_components];
    for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
      alpha_e[i] = 1.;
    // Init u_i_min, u_i_max
    double u_i_min[Equations<equationsType, dim>::n_components];
    double u_i_max[Equations<equationsType, dim>::n_components];
    for (int k = 0; k < Equations<equationsType, dim>::n_components; k++)
    {
      u_i_min[k] = u_c[k];
      u_i_max[k] = u_c[k];
    }

    for (unsigned int vertex_i = 0; vertex_i < GeometryInfo<dim>::vertices_per_cell; ++vertex_i)
    {
      // For all vertices -> v_i
      for (auto dof_indices_neighbor : data->neighbor_dof_indices[vertex_i])
      {
        if (dof_indices_neighbor.size() == 0)
          continue;
        bool u_i_extrema_set[Equations<equationsType, dim>::n_components];
        for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
          u_i_extrema_set[i] = false;

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          if (is_primitive[i])
          {
            // Here we rely on the fact, that the constant basis fn is the first one.
            if (!u_i_extrema_set[component_ii[i]])
            {
              double val = current_unlimited_solution(dof_indices_neighbor[i]);
              if (parameters.debug)
              {
                if (val < u_i_min[component_ii[i]])
                  std::cout << "\tdecreasing u_i_min to: " << val << std::endl;
                if (val > u_i_max[component_ii[i]])
                  std::cout << "\tincreasing u_i_max to: " << val << std::endl;
              }
              u_i_min[component_ii[i]] = std::min(u_i_min[component_ii[i]], val);
              u_i_max[component_ii[i]] = std::max(u_i_max[component_ii[i]], val);
              u_i_extrema_set[component_ii[i]] = true;
            }
          }
        }
      }
    }

    // Based on u_i_min, u_i_max, u_i, get alpha_e
    for (unsigned int vertex_i = 0; vertex_i < GeometryInfo<dim>::vertices_per_cell; ++vertex_i)
    {
      // (!!!) Find out u_i
      Vector<double> u_i(Equations<equationsType, dim>::n_components);
      const Point<dim> p_cell = mapping.transform_real_to_unit_cell(cell, data->vertexPoint[vertex_i]);
      const Quadrature<dim> one_point_quadrature(GeometryInfo<dim>::project_to_unit_cell(p_cell));
      FEValues<dim> fe_values(mapping, fe, one_point_quadrature, update_values);
      fe_values.reinit(cell);
      std::vector<Vector<double> > u_value(1, Vector<double>(fe.n_components()));
      fe_values.get_function_values(current_unlimited_solution, u_value);
      u_i = u_value[0];

      if (parameters.debug)
      {
        std::cout << "\tv_i: " << cell->vertex(vertex_i) << ", values: ";
        for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
          std::cout << u_i[i] << (i == Equations<equationsType, dim>::n_components - 1 ? "" : ", ");
        std::cout << std::endl;
      }

      for (int k = 0; k < Equations<equationsType, dim>::n_components; k++)
        if (std::abs((u_c[k] - u_i[k]) / u_c[k]) > NEGLIGIBLE)
        {
          alpha_e[k] = std::min(alpha_e[k], ((u_i[k] - u_c[k]) > 0.) ? std::min(1.0, (u_i_max[k] - u_c[k]) / (u_i[k] - u_c[k])) : std::min(1.0, (u_i_min[k] - u_c[k]) / (u_i[k] - u_c[k])));
          if (parameters.debug)
            std::cout << "\talpha_e[" << k << "]: " << alpha_e[k] << std::endl;
        }
    }

    for (int k = 0; k < Equations<equationsType, dim>::n_components; k++)
      for (int i = 0; i < data->lambda_indices_to_multiply[k].size(); i++)
        current_limited_solution(data->lambda_indices_to_multiply[k][i]) *= alpha_e[k];

    double alpha_e_B = std::min(std::min(alpha_e[5], alpha_e[6]), alpha_e[7]);
    for (int i = 0; i < data->lambda_indices_to_multiply_all_B_components.size(); i++)
      current_limited_solution(data->lambda_indices_to_multiply_all_B_components[i]) *= alpha_e_B;
  }
}

template class SlopeLimiter<EquationsTypeMhd, 3>;
template class VertexBasedSlopeLimiter<EquationsTypeMhd, 3>;
template class BarthJespersenSlopeLimiter<EquationsTypeMhd, 3>;
