/*
 * DealiiExtensions.cpp
 *
 *  Created on: Dec 11, 2014
 *      Author: kraemer
 */

#include "DealiiExtensions.h"

#include <deal.II/base/thread_management.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <numeric>

using namespace dealii;
using namespace DoFTools;

namespace DealIIExtensions
{

  template<class DH, class SparsityPattern>
  void make_sparser_flux_sparsity_pattern(const DH &dof, SparsityPattern &sparsity, const ConstraintMatrix &constraints, const std::vector<std::array<int, 3> >& boundaries,
    const PeriodicCellMap<DH::dimension>& cell_map, FEFaceValues<DH::dimension>* fe_face, const bool keep_constrained_dofs, const types::subdomain_id subdomain_id)
  {
#ifdef DEBUG
    const types::global_dof_index n_dofs = dof.n_dofs();
    AssertDimension(sparsity.n_rows(), n_dofs);
    AssertDimension(sparsity.n_cols(), n_dofs);
#endif

    if (fe_face != NULL)
    {
      Assert(fe_face->get_fe().has_support_points(),
        ExcMessage("Sparser flux sparsity pattern makes only sense for elements with support points"));
    }

    // If we have a distributed::Mesh only allow locally_owned
    // subdomain. Not setting a subdomain is also okay, because we skip
    // ghost cells in the loop below.
    Assert(
      (dof.get_tria().locally_owned_subdomain() == numbers::invalid_subdomain_id) || (subdomain_id == numbers::invalid_subdomain_id) || (subdomain_id == dof.get_tria().locally_owned_subdomain()),
      ExcMessage("For parallel::distributed::Mesh objects and " "associated DoF handler objects, asking for any subdomain other " "than the locally owned one does not make sense."));

    std::vector<types::global_dof_index> dofs_on_this_cell;
    std::vector<types::global_dof_index> dofs_on_other_cell;
    std::vector<types::global_dof_index> dofs_on_this_face;
    std::vector<types::global_dof_index> dofs_on_other_face;
    dofs_on_this_cell.reserve(max_dofs_per_cell(dof));
    dofs_on_other_cell.reserve(max_dofs_per_cell(dof));
    dofs_on_this_face.reserve(max_dofs_per_cell(dof));
    dofs_on_other_face.reserve(max_dofs_per_cell(dof));

    // TODO: in an old implementation, we used user flags before to tag
    // faces that were already touched. this way, we could reduce the work
    // a little bit. now, we instead add only data from one side. this
    // should be OK, but we need to actually verify it.

    // In case we work with a distributed sparsity pattern of Trilinos
    // type, we only have to do the work if the current cell is owned by
    // the calling processor. Otherwise, just continue.
    for (typename DH::active_cell_iterator cell = dof.begin_active(); cell != dof.end(); ++cell)
    {
      if (((subdomain_id == numbers::invalid_subdomain_id) || (subdomain_id == cell->subdomain_id())) && cell->is_locally_owned())
      {
        const unsigned int n_dofs_on_this_cell = cell->get_fe().dofs_per_cell;
        dofs_on_this_cell.resize(n_dofs_on_this_cell);
        cell->get_dof_indices(dofs_on_this_cell);

        // make sparsity pattern for this cell. if no constraints pattern
        // was given, then the following call acts as if simply no
        // constraints existed
        constraints.add_entries_local_to_global(dofs_on_this_cell, sparsity, keep_constrained_dofs);

        for (unsigned int face = 0; face < GeometryInfo<DH::dimension>::faces_per_cell; ++face)
        {
          typename DH::face_iterator this_face = cell->face(face);
          typename DH::face_iterator other_face;
          typename DH::cell_iterator neighbor(cell);
          unsigned int neighbor_face = 1000;
          if (cell->at_boundary(face))
          {
            for (int p = 0; p < boundaries.size(); p++)
            {
              if (boundaries[p][0] == this_face->boundary_id() || boundaries[p][1] == this_face->boundary_id())
              {
                const FacePair<DH::dimension>& face_pair = cell_map.find(cell)->second;
                if (cell_map.find(cell) == cell_map.end())
                {
                  std::cout << "Something wrong (unable to find in cell_map)" << std::endl;
                  continue;
                }
                neighbor = ((*(face_pair.cell[0])).active_cell_index() == (*cell).active_cell_index()) ? face_pair.cell[1] : face_pair.cell[0];
                neighbor_face = ((*(face_pair.cell[0])).active_cell_index() == (*cell).active_cell_index()) ? face_pair.face_idx[1] : face_pair.face_idx[0];
                other_face = neighbor->face(neighbor_face);
                break;
              }
            }
            if (neighbor_face == 1000)
            {
              continue;
            }
          }
          else {
            neighbor = cell->neighbor(face);
            neighbor_face = 1000; // indicator that it has to be assembled, later
          }

          // specify whether pairwise coupling is valid
          bool pairwise_coupling_valid = (fe_face != NULL);

          dofs_on_this_face.clear();
          dofs_on_other_face.clear();
          if (!pairwise_coupling_valid)
          {
            // get face dofs
            for (size_t i = 0; i < n_dofs_on_this_cell; i++)
            {
              if (cell->get_fe().has_support_on_face(i, face))
              {
                dofs_on_this_face.push_back(
                  dofs_on_this_cell.at(i));
              }
            }
          }

          if (neighbor->has_children())
          {
            for (unsigned int sub_nr = 0;
              sub_nr != this_face->number_of_children();
              ++sub_nr)
            {
              const typename DH::cell_iterator sub_neighbor =
                cell->neighbor_child_on_subface(face, sub_nr);

              const unsigned int n_dofs_on_neighbor =
                sub_neighbor->get_fe().dofs_per_cell;
              dofs_on_other_cell.resize(n_dofs_on_neighbor);
              sub_neighbor->get_dof_indices(dofs_on_other_cell);

              // identify which sub_neighbor face is child to this face
              unsigned int sub_neighbor_face;
              for (sub_neighbor_face = 0;
                sub_neighbor_face
                < GeometryInfo<DH::dimension>::faces_per_cell;
                ++sub_neighbor_face)
              {
                other_face = sub_neighbor->face(sub_neighbor_face);
                if (sub_neighbor->neighbor(sub_neighbor_face)
                  == cell)
                {
                  break;
                }
                Assert(
                  sub_neighbor_face + 1 < GeometryInfo<DH::dimension>::faces_per_cell,
                  ExcMessage("Neighbor face was not found, but needed for constructing the sparsity pattern"));
              }

              // Couple all dofs on common face
              dofs_on_other_face.clear();
              for (size_t i = 0; i < n_dofs_on_neighbor; i++)
              {
                if (sub_neighbor->get_fe().has_support_on_face(i,
                  sub_neighbor_face))
                {
                  dofs_on_other_face.push_back(
                    dofs_on_other_cell.at(i));
                }
              }
              Assert(
                dofs_on_this_face.size() * dofs_on_other_face.size() > 0,
                ExcMessage("Size of at least one dof vector is 0."));

              // Add entries to sparsity pattern
              constraints.add_entries_local_to_global(
                dofs_on_this_face, dofs_on_other_face, sparsity,
                keep_constrained_dofs);
              constraints.add_entries_local_to_global(
                dofs_on_other_face, dofs_on_this_face, sparsity,
                keep_constrained_dofs);
              // only need to add this when the neighbor is not
              // owned by the current processor, otherwise we add
              // the entries for the neighbor there
              if (sub_neighbor->subdomain_id()
                != cell->subdomain_id())
                constraints.add_entries_local_to_global(
                  dofs_on_other_cell, sparsity,
                  keep_constrained_dofs);
            }
          }
          else {
            // Refinement edges are taken care of by coarser
            // cells
            if (!cell->at_boundary(face))
            {
              if (cell->neighbor_is_coarser(face)
                && neighbor->subdomain_id()
                == cell->subdomain_id())
                continue;
            }
            const unsigned int n_dofs_on_neighbor =
              neighbor->get_fe().dofs_per_cell;
            dofs_on_other_cell.resize(n_dofs_on_neighbor);
            neighbor->get_dof_indices(dofs_on_other_cell);

            // identify which neighbor face belongs to this face
            if (neighbor_face == 1000)
            {
              for (neighbor_face = 0; neighbor_face < GeometryInfo<DH::dimension>::faces_per_cell; ++neighbor_face)
              {
                other_face = neighbor->face(neighbor_face);
                if (*other_face == *this_face)
                  break;
                Assert(neighbor_face + 1 < GeometryInfo<DH::dimension>::faces_per_cell, ExcMessage("Neighbor face was not found, but needed for constructing the sparsity pattern"));
              }
            }

            if (!pairwise_coupling_valid)
            {
              // Method 1) Couple all dofs on common face
              dofs_on_other_face.clear();
              for (size_t i = 0; i < n_dofs_on_neighbor; i++)
              {
                if (neighbor->get_fe().has_support_on_face(i, neighbor_face))
                  dofs_on_other_face.push_back(dofs_on_other_cell.at(i));
              }
              Assert(
                dofs_on_this_face.size() * dofs_on_other_face.size() > 0, ExcMessage("Size of at least one dof vector is 0."));

              // Add entries to sparsity pattern
              constraints.add_entries_local_to_global(dofs_on_this_face, dofs_on_other_face, sparsity, keep_constrained_dofs);
            }
            else {
              // Method 2) if possible: make unique relation between neighboring dofs
              fe_face->reinit(cell, face);
              // bring dofs at this face into right order
              for (size_t i = 0; i < n_dofs_on_this_cell; i++)
              {
                int unique = 0;
                for (size_t q = 0; q < fe_face->n_quadrature_points; q++)
                {
                  if (fe_face->shape_value(i, q) > 1e-10)
                  {
                    unique += 1;
                    dofs_on_this_face.push_back(dofs_on_this_cell.at(i));
                  }
                }
                // Test, if the relationship doF <-> quadrature points in unique
                assert(unique <= 1);
              }
              // bring dofs at other face in right order
              fe_face->reinit(neighbor, neighbor_face);
              for (size_t i = 0; i < n_dofs_on_neighbor; i++)
              {
                int unique = 0;
                for (size_t q = 0; q < fe_face->n_quadrature_points; q++)
                {
                  if (fe_face->shape_value(i, q) > 1e-10)
                  {
                    unique += 1;
                    dofs_on_other_face.push_back(dofs_on_other_cell.at(i));
                  }
                }
                // Test, if the relationship doF <-> quadrature points in unique
                assert(unique <= 1);
              }

              if (2 == DH::dimension)
              {
                AssertDimension(dofs_on_this_face.size(), sqrt(n_dofs_on_this_cell));
              }
              else if (3 == DH::dimension)
              {
                AssertDimension(dofs_on_this_face.size()* sqrt(dofs_on_this_face.size()), n_dofs_on_this_cell);
              }
              AssertDimension(dofs_on_other_face.size(), dofs_on_this_face.size());

              // couple only individual dofs with each other
              std::vector<types::global_dof_index> dof_this(1);
              std::vector<types::global_dof_index> dof_other(1);
              for (size_t i = 0; i < dofs_on_this_face.size(); i++)
              {
                dof_this.at(0) = dofs_on_this_face.at(i);
                dof_other.at(0) = dofs_on_other_face.at(i);

                // Add entries to sparsity pattern
                constraints.add_entries_local_to_global(dof_this, dof_other, sparsity, keep_constrained_dofs);
              }
            }
            // only need to add these in case the neighbor cell
            // is not locally owned - otherwise, we touch each
            // face twice and hence put the indices the other way
            // around
            if (!(neighbor->active())
              || (neighbor->subdomain_id() != cell->subdomain_id()))
            {
              if (!pairwise_coupling_valid)
              {
                constraints.add_entries_local_to_global(dofs_on_other_face, dofs_on_this_face, sparsity, keep_constrained_dofs);
              }
              else
              {
                // couple only individual dofs with each other
                std::vector<types::global_dof_index> dof_this(1);
                std::vector<types::global_dof_index> dof_other(1);
                for (size_t i = 0; i < dofs_on_this_face.size(); i++)
                {
                  dof_this.at(0) = dofs_on_this_face.at(i);
                  dof_other.at(0) = dofs_on_other_face.at(i);

                  // Add entries to sparsity pattern
                  constraints.add_entries_local_to_global(dof_this, dof_other, sparsity, keep_constrained_dofs);
                }
              }
              if (neighbor->subdomain_id() != cell->subdomain_id())
                constraints.add_entries_local_to_global(dofs_on_other_cell, sparsity, keep_constrained_dofs);
            }
          }
        }
      }
    }
  }

  template<typename DH>
  void make_periodicity_map_dg(const typename DH::cell_iterator &cell_1, const typename identity<typename DH::cell_iterator>::type &cell_2, size_t face_nr_1, size_t face_nr_2,
    PeriodicCellMap<DH::dimension>& cell_map, const bool face_orientation, const bool face_flip, const bool face_rotation)
  {
    static const int dim = DH::dimension;
    static const int spacedim = DH::space_dimension;

    typedef typename DH::face_iterator FaceIterator;
    FaceIterator face_1 = cell_1->face(face_nr_1);
    FaceIterator face_2 = cell_2->face(face_nr_2);

    Assert(
      (dim != 1) || (face_orientation == true && face_flip == false && face_rotation == false),
      ExcMessage("The supplied orientation " "(face_orientation, face_flip, face_rotation) " "is invalid for 1D"));

    Assert((dim != 2) || (face_orientation == true && face_rotation == false),
      ExcMessage("The supplied orientation " "(face_orientation, face_flip, face_rotation) " "is invalid for 2D"));

    Assert(face_1 != face_2,
      ExcMessage("face_1 and face_2 are equal! Cannot create boundary conditions DoFs " "on the very same face"));

    Assert(face_1->at_boundary() && face_2->at_boundary(),
      ExcMessage("Faces for periodicity constraints must be on the " "boundary"));

    // A lookup table on how to go through the child cells depending on the
    // orientation:
    // see Documentation of GeometryInfo for details

    static const int lookup_table_2d[2][2] =
    {
      //          flip:
          { 0, 1 }, //  false
          { 1, 0 }, //  true
    };

    static const int lookup_table_3d[2][2][2][4] =
    {
      //                    orientation flip  rotation
          { { { 0, 2, 1, 3 }, //  false       false false
              { 2, 3, 0, 1 }, //  false       false true
              }, { { 3, 1, 2, 0 }, //  false       true  false
                  { 1, 0, 3, 2 }, //  false       true  true
              }, }, { { { 0, 1, 2, 3 }, //  true        false false
              { 1, 3, 0, 2 }, //  true        false true
              }, { { 3, 2, 1, 0 }, //  true        true  false
                  { 2, 0, 3, 1 }, //  true        true  true
              }, }, };

    if (cell_1->has_children() && cell_2->has_children())
    {
      // In the case that both faces have children, we loop over all
      // children and apply make_periodicty_constrains recursively:
      Assert(
        face_1->n_children() == GeometryInfo<dim>::max_children_per_face && face_2->n_children() == GeometryInfo<dim>::max_children_per_face,
        ExcNotImplemented());

      for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_face; ++i)
      {
        // Lookup the index for the second face
        unsigned int j;
        switch (dim)
        {
        case 2:
          j = lookup_table_2d[face_flip][i];
          break;
        case 3:
          j =
            lookup_table_3d[face_orientation][face_flip][face_rotation][i];
          break;
        default:
          AssertThrow(false, ExcNotImplemented())
            ;
        }

        // find subcell ids that belong to the subface indices
        size_t child_cell_1 = GeometryInfo<dim>::child_cell_on_face(
          cell_1->refinement_case(), face_nr_1, i, face_orientation,
          face_flip, face_rotation, face_1->refinement_case());
        size_t child_cell_2 = GeometryInfo<dim>::child_cell_on_face(
          cell_2->refinement_case(), face_nr_2, j, face_orientation,
          face_flip, face_rotation, face_2->refinement_case());

        // precondition: subcell has the same orientation as cell (so that the face numbers coincide)
        // recursive call
        make_periodicity_map_dg<DH>(cell_1->child(child_cell_1),
          cell_2->child(child_cell_2), face_nr_1, face_nr_2, cell_map,
          face_orientation, face_flip, face_rotation);
      }
    }
    else
    {
      // Otherwise at least one of the two faces is active and
      // we need to do some work and enter the periodic face pairs!

      // This case could only be allowed with anisotropic refinement,
      // because the coarser cell is not allowed to have subfaces at the boundary.
      // Anisotropic refinement is also left out here for simplicity.
      // Consequently, opposite cells at periodic boundaries have to have the same
      // refinement level.
      if (cell_1->has_children())
        Assert((cell_1->has_children() && cell_2->has_children()) || (!cell_1->has_children() && !cell_2->has_children()),
          ExcMessage("Refinement levels of the matched cells must match."));

      // taken from make_periodicity_constraints: make sure faces are not artificial
      const unsigned int face_1_index = face_1->nth_active_fe_index(0);
      const unsigned int face_2_index = face_2->nth_active_fe_index(0);
      const unsigned int dofs_per_face = face_1->get_fe(face_1_index).dofs_per_face;
      std::vector<types::global_dof_index> dofs_1(dofs_per_face);
      std::vector<types::global_dof_index> dofs_2(dofs_per_face);
      face_1->get_dof_indices(dofs_1, face_1_index);
      face_2->get_dof_indices(dofs_2, face_2_index);
      for (unsigned int i = 0; i < dofs_per_face; i++)
      {
        if (dofs_1[i] == numbers::invalid_dof_index
          || dofs_2[i] == numbers::invalid_dof_index)
        {
          /* If either of these faces have no indices, stop.  This is so
           * that there is no attempt to match artificial cells of
           * parallel distributed triangulations.
           *
           * While it seems like we ought to be able to avoid even calling
           * set_periodicity_constraints for artificial faces, this
           * situation can arise when a face that is being made periodic
           * is only partially touched by the local subdomain.
           * make_periodicity_constraints will be called recursively even
           * for the section of the face that is not touched by the local
           * subdomain.
           *
           * Until there is a better way to determine if the cells that
           * neighbor a face are artificial, we simply test to see if the
           * face does not have a valid dof initialization.
           */
          return;
        }
      }

      // make sure cells are at least ghost cells (or even locally owned)
      if ((cell_1->has_children()) || (cell_2->has_children()))
        return;
      if ((cell_1->is_artificial()) || (cell_2->is_artificial()))
        return;
      /*Assert(not cell_1->is_artificial(),
          ExcMessage ("Cell at periodic boundary must not be artificial."));
      Assert(not cell_2->is_artificial(),
          ExcMessage ("Cell at periodic boundary must not be artificial."));
       */

       // insert periodic face pair for both cells
      dealii::GridTools::PeriodicFacePair<dealii::TriaIterator<dealii::CellAccessor<dim, spacedim> > > face_pair;
      face_pair.cell[0] = cell_1;
      face_pair.cell[1] = cell_2;
      face_pair.face_idx[0] = face_nr_1;
      face_pair.face_idx[1] = face_nr_2;
      face_pair.orientation[0] = face_orientation;
      face_pair.orientation[1] = face_flip;
      face_pair.orientation[2] = face_rotation;
      cell_map.insert(
        std::pair<
        dealii::TriaIterator<dealii::CellAccessor<dim, spacedim> >,
        dealii::GridTools::PeriodicFacePair<
        dealii::TriaIterator<
        dealii::CellAccessor<dim, spacedim> > > >(
          cell_1, face_pair));

      cell_map.insert(
        std::pair<
        dealii::TriaIterator<dealii::CellAccessor<dim, spacedim> >,
        dealii::GridTools::PeriodicFacePair<
        dealii::TriaIterator<
        dealii::CellAccessor<dim, spacedim> > > >(
          cell_2, face_pair));
    }
  }

  template<typename DH>
  void make_periodicity_map_dg(
    const std::vector<
    GridTools::PeriodicFacePair<typename DH::cell_iterator> > &periodic_faces,
    PeriodicCellMap<DH::dimension>& cell_map)
  {
    typedef std::vector<GridTools::PeriodicFacePair<typename DH::cell_iterator> > FaceVector;
    typename FaceVector::const_iterator it, end_periodic;
    it = periodic_faces.begin();
    end_periodic = periodic_faces.end();

    // Loop over all periodic faces...
    for (; it != end_periodic; ++it)
    {
      const typename DH::cell_iterator cell_1 = it->cell[0];
      size_t face_id_1 = it->face_idx[0];
      const typename DH::cell_iterator cell_2 = it->cell[1];
      size_t face_id_2 = it->face_idx[1];

      Assert(
        cell_1->face(face_id_1)->at_boundary() && cell_2->face(face_id_2)->at_boundary(),
        ExcInternalError());

      Assert(cell_1->face(face_id_1) != cell_2->face(face_id_2),
        ExcInternalError());

      // ... and apply the low level function to
      // every matching pair:
      make_periodicity_map_dg<DH>(cell_1, cell_2, face_id_1, face_id_2,
        cell_map, it->orientation[0], it->orientation[1],
        it->orientation[2]);
    }
  }

  // High level interface variants:

  template<typename DH>
  void make_periodicity_map_dg(const DH &dof_handler, size_t b_id1, size_t b_id2,
    const int direction, PeriodicCellMap<DH::dimension>& cell_map)
  {
    static const int space_dim = DH::space_dimension;
    (void)space_dim;
    Assert(0 <= direction && direction < space_dim,
      ExcIndexRange(direction, 0, space_dim));

    Assert(b_id1 != b_id2,
      ExcMessage("The boundary indicators b_id1 and b_id2 must be" "different to denote different boundaries."));

    std::vector<GridTools::PeriodicFacePair<typename DH::cell_iterator> > matched_pairs;

    // Collect matching periodic cells on the coarsest level:
    GridTools::collect_periodic_faces(dof_handler, b_id1, b_id2, direction, matched_pairs);

    // call function to make the periodicity map
    make_periodicity_map_dg<DH>(matched_pairs, cell_map);
  }

  template<class DH>
  void extract_dofs_with_support_on_boundary(const DH &dof_handler,
    const ComponentMask &component_mask, std::vector<bool> &selected_dofs,
    const std::set<types::boundary_id> &boundary_ids)
  {
    Assert(component_mask.represents_n_components(n_components(dof_handler)),
      ExcMessage("This component mask has the wrong size."));
    Assert(
      boundary_ids.find(numbers::internal_face_boundary_id) == boundary_ids.end(),
      ExcInvalidBoundaryIndicator());
    // let's see whether we have to check for certain boundary indicators
    // || whether we can accept all
    const bool check_boundary_id = (boundary_ids.size() != 0);

    // also see whether we have to check whether a certain vector component
    // is selected, || all
    const bool check_vector_component =
      (component_mask.represents_the_all_selected_mask() == false);

    // clear and reset array by default values
    selected_dofs.clear();
    selected_dofs.resize(dof_handler.n_dofs(), false);
    std::vector<types::global_dof_index> cell_dof_indices;
    cell_dof_indices.reserve(max_dofs_per_cell(dof_handler));

    // now loop over all cells and check whether their faces are at the
    // boundary. note that we need not take special care of single lines
    // being at the boundary (using @p{cell->has_boundary_lines}), since we
    // do not support boundaries of dimension dim-2, and so every isolated
    // boundary line is also part of a boundary face which we will be
    // visiting sooner || later
    for (typename DH::active_cell_iterator cell = dof_handler.begin_active();
      cell != dof_handler.end(); ++cell)
    {
      if (cell->is_locally_owned())
      {
        for (unsigned int face = 0;
          face < GeometryInfo<DH::dimension>::faces_per_cell;
          ++face)
        {
          if (cell->at_boundary(face))
          {
            if (!check_boundary_id
              || (boundary_ids.find(
                cell->face(face)->boundary_id())
                != boundary_ids.end()))
            {
              const FiniteElement<DH::dimension, DH::space_dimension> &fe =
                cell->get_fe();

              const unsigned int dofs_per_cell = fe.dofs_per_cell;
              cell_dof_indices.resize(dofs_per_cell);
              cell->get_dof_indices(cell_dof_indices);

              for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
              {
                if (fe.has_support_on_face(i, face))
                {
                  if (!check_vector_component)
                    selected_dofs[cell_dof_indices[i]] = true;
                  else
                    // check for component is required. somewhat tricky
                    // as usual for the case that the shape function is
                    // non-primitive, but use usual convention (see docs)
                  {
                    if (fe.is_primitive(i))
                      selected_dofs[cell_dof_indices[i]] =
                      (component_mask[fe.system_to_component_index(
                        i).first] == true);
                    else // not primitive
                    {
                      const unsigned int first_nonzero_comp =
                        fe.get_nonzero_components(i).first_selected_component();
                      Assert(
                        first_nonzero_comp < fe.n_components(),
                        ExcInternalError());
                      selected_dofs[cell_dof_indices[i]] =
                        (component_mask[first_nonzero_comp]
                          == true);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  typedef TrilinosWrappers::SparsityPattern SP;
  //for (SP : SPARSITY_PATTERNS; deal_II_dimension : DIMENSIONS)
  //for (size_t deal_II_dimension = 1; deal_II_dimension < 4; deal_II_dimension++)
  //{
  template void
    make_sparser_flux_sparsity_pattern<DoFHandler<2>, SP>(const DoFHandler<2> &dof,
      SP &sparsity, const ConstraintMatrix &constraints,
      const std::vector<std::array<int, 3> >& boundaries,
      const PeriodicCellMap<2>& cell_map, FEFaceValues<2>* fe_face,
      const bool, const unsigned int);

  template void
    make_sparser_flux_sparsity_pattern<DoFHandler<3>, SP>(const DoFHandler<3> &dof,
      SP &sparsity, const ConstraintMatrix &constraints,
      const std::vector<std::array<int, 3> >& boundaries,
      const PeriodicCellMap<3>& cell_map, FEFaceValues<3>* fe_face,
      const bool, const unsigned int);

  template void
    make_sparser_flux_sparsity_pattern<DoFHandler<3>, dealii::DynamicSparsityPattern>(const DoFHandler<3> &dof,
      dealii::DynamicSparsityPattern &sparsity, const ConstraintMatrix &constraints,
      const std::vector<std::array<int, 3> >& boundaries,
      const PeriodicCellMap<3>& cell_map, FEFaceValues<3>* fe_face,
      const bool, const unsigned int);

  template
    void make_periodicity_map_dg<DoFHandler<2> >(
      const typename DoFHandler<2>::cell_iterator &cell_1,
      const identity<typename DoFHandler<2>::cell_iterator>::type &cell_2,
      size_t face_nr_1, size_t face_nr_2, PeriodicCellMap<2>& cell_map,
      const bool face_orientation, const bool face_flip,
      const bool face_rotation);
  template
    void make_periodicity_map_dg<DoFHandler<3> >(
      const typename DoFHandler<3>::cell_iterator &cell_1,
      const identity<typename DoFHandler<3>::cell_iterator>::type &cell_2,
      size_t face_nr_1, size_t face_nr_2, PeriodicCellMap<3>& cell_map,
      const bool face_orientation, const bool face_flip,
      const bool face_rotation);

  template
    void make_periodicity_map_dg<DoFHandler<2> >(
      const std::vector<
      GridTools::PeriodicFacePair<
      typename DoFHandler<2>::cell_iterator> > &periodic_faces,
      PeriodicCellMap<2>& cell_map);

  template
    void make_periodicity_map_dg<DoFHandler<3> >(
      const std::vector<
      GridTools::PeriodicFacePair<
      typename DoFHandler<3>::cell_iterator> > &periodic_faces,
      PeriodicCellMap<3>& cell_map);

  template
    void make_periodicity_map_dg<DoFHandler<2> >(const DoFHandler<2> &dof_handler,
      size_t b_id1, size_t b_id2, const int direction,
      PeriodicCellMap<2>& cell_map);
  template
    void make_periodicity_map_dg<DoFHandler<3> >(const DoFHandler<3> &dof_handler,
      size_t b_id1, size_t b_id2, const int direction,
      PeriodicCellMap<3>& cell_map);
  //}s

  template
    void extract_dofs_with_support_on_boundary(
      const dealii::DoFHandler<2> &dof_handler,
      const ComponentMask &component_mask, std::vector<bool> &selected_dofs,
      const std::set<types::boundary_id> &boundary_ids);

  template
    void extract_dofs_with_support_on_boundary(
      const dealii::DoFHandler<3> &dof_handler,
      const ComponentMask &component_mask, std::vector<bool> &selected_dofs,
      const std::set<types::boundary_id> &boundary_ids);

}
