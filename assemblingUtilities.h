#ifndef _ASSEMBLING_UTILITIES_H
#define _ASSEMBLING_UTILITIES_H

#include "problem.h"

template <EquationsType equationsType, int dim> class Problem;

// Class that accepts all input from the user, provides interface for output, etc.
// Should not be changed.
template <EquationsType equationsType, int dim>
class AssemblingUtilities
{
public:
  void set_problem(Problem<equationsType, dim>* problem);
  bool AssemblingUtilities<equationsType, dim>::is_refinement_compatible_with_subface(int face_no, RefinementCase<dim> ref_case, int child, int subface_no, RefinementCase<dim - 1> face_ref_case);
  bool AssemblingUtilities<equationsType, dim>::is_refinement_within_current_cell(int face_no, RefinementCase<dim> ref_case, int child);

private:
  Problem<equationsType, dim>* problem;
};
#endif