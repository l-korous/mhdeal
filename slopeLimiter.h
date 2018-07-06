#ifndef _SLOPE_LIMITER_H
#define _SLOPE_LIMITER_H

#include "equationsMhd.h"
#include "parameters.h"

template <EquationsType equationsType, int dim>
class SlopeLimiter
{
public:
  SlopeLimiter(const Parameters<dim>& parameters, const MappingQ1<dim>& mapping, const FESystem<dim>& fe, DoFHandler<dim>& dof_handler, 
     unsigned int& dofs_per_cell,
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<dim>& triangulation,
#else
    Triangulation<dim>& triangulation,
#endif
    std::vector<types::global_dof_index>& dof_indices, std::array <unsigned short, BASIS_FN_COUNT>& component_ii, std::array <bool, BASIS_FN_COUNT>& is_primitive) : 
    parameters(parameters),
    mapping(mapping),
    fe(fe),
    dof_handler(dof_handler),
    dofs_per_cell(dofs_per_cell),
    triangulation(triangulation),
    dof_indices(dof_indices),
    component_ii(component_ii),
    is_primitive(is_primitive)
    {};

  // Not const because of caching.
  virtual void postprocess(TrilinosWrappers::MPI::Vector& current_limited_solution, TrilinosWrappers::MPI::Vector& current_unlimited_solution) = 0;
  virtual void flush_cache() = 0;
protected:
  // Performs a single global assembly.
  struct PostprocessData
  {
    Point<dim> center;
    unsigned int vertexIndex[GeometryInfo<dim>::vertices_per_cell];
    Point<dim> vertexPoint[GeometryInfo<dim>::vertices_per_cell];
    std::vector<unsigned int> lambda_indices_to_multiply[Equations<equationsType, dim>::n_components];
    std::vector<unsigned int> lambda_indices_to_multiply_all_B_components;
    // TODO This may not work for adaptivity (max. number of cells sharing a vertes might be higher than 8)
    std::vector<types::global_dof_index> neighbor_dof_indices[GeometryInfo<dim>::vertices_per_cell][GeometryInfo<dim>::vertices_per_cell];
    unsigned short neighbor_count;
    std::array<bool, GeometryInfo<dim>::vertices_per_cell> vertex_is_at_nonperiodic_boundary;
  };

  std::map<unsigned int, PostprocessData> postprocessData;
  
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& triangulation;
#else
  Triangulation<dim>& triangulation;
#endif
  const Parameters<dim>& parameters;
  const MappingQ1<dim>& mapping;
  const FESystem<dim>& fe;
  DoFHandler<dim>& dof_handler;
  unsigned int& dofs_per_cell;
  std::vector<types::global_dof_index>& dof_indices;
  std::array <unsigned short, BASIS_FN_COUNT>& component_ii;
  std::array <bool, BASIS_FN_COUNT>& is_primitive;
};

template <EquationsType equationsType, int dim>
class VertexBasedSlopeLimiter : public SlopeLimiter<equationsType, dim>
{
public:
  VertexBasedSlopeLimiter(const Parameters<dim>& parameters, const MappingQ1<dim>& mapping, const FESystem<dim>& fe, DoFHandler<dim>& dof_handler,  
    unsigned int& dofs_per_cell,
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<dim>& triangulation,
#else
    Triangulation<dim>& triangulation,
#endif
    std::vector<types::global_dof_index>& dof_indices, std::array <unsigned short, BASIS_FN_COUNT>& component_ii, std::array <bool, BASIS_FN_COUNT>& is_primitive) : 
    SlopeLimiter<equationsType, dim>(parameters, mapping, fe, dof_handler, dofs_per_cell, triangulation, dof_indices, component_ii, is_primitive) {};
  // Not const because of caching.
  virtual void postprocess(TrilinosWrappers::MPI::Vector& current_limited_solution, TrilinosWrappers::MPI::Vector& current_unlimited_solution);
  void flush_cache();
};

template <EquationsType equationsType, int dim>
class BarthJespersenSlopeLimiter : public SlopeLimiter<equationsType, dim>
{
public:
  BarthJespersenSlopeLimiter(const Parameters<dim>& parameters, const MappingQ1<dim>& mapping, const FESystem<dim>& fe, DoFHandler<dim>& dof_handler,  
    unsigned int& dofs_per_cell,
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<dim>& triangulation,
#else
    Triangulation<dim>& triangulation,
#endif
    std::vector<types::global_dof_index>& dof_indices, std::array <unsigned short, BASIS_FN_COUNT>& component_ii, std::array <bool, BASIS_FN_COUNT>& is_primitive) : 
    SlopeLimiter<equationsType, dim>(parameters, mapping, fe, dof_handler, dofs_per_cell, triangulation, dof_indices, component_ii, is_primitive) {};
  // Not const because of caching.
  virtual void postprocess(TrilinosWrappers::MPI::Vector& current_limited_solution, TrilinosWrappers::MPI::Vector& current_unlimited_solution);
  void flush_cache();
};

#endif
