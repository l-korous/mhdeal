#include "assemblingUtilities.h"

template <EquationsType equationsType, int dim>
void AssemblingUtilities<equationsType, dim>::set_problem(Problem<equationsType, dim>* problem_)
{
  this->problem = problem_;
}

/*
* @verbatim
*       *-------*        *-------*
*      /|       |       /       /|
*     / |   3   |      /   5   / |
*    /  |       |     /       /  |
*   *   |       |    *-------*   |
*   | 0 *-------*    |       | 1 *
*   |  /       /     |       |  /
*   | /   4   /      |   2   | /
*   |/       /       |       |/
*   *-------*        *-------*
* @endverbatim

* @verbatim
*            RefinementCase<3>::cut_xy
*
*       *----*----*        *----*----*
*      /|    |    |       / 2  /  3 /|
*     * |    |    |      *----*----* |
*    /| | 2  |  3 |     / 0  /  1 /| |
*   * |2|    |    |    *----*----* |3|
*   | | |    |    |    |    |    | | |
*   |0| *----*----*    |    |    |1| *
*   | |/ 2  /  3 /     | 0  |  1 | |/
*   | *----*----*      |    |    | *
*   |/ 0  /  1 /       |    |    |/
*   *----*----*        *----*----*
* @endverbatim
*
* @verbatim
*            RefinementCase<3>::cut_xyz
*       *----*----*        *----*----*
*      /| 6  |  7 |       / 6  /  7 /|
*     *6|    |    |      *----*----*7|
*    /| *----*----*     / 4  /  5 /| *
*   * |/|    |    |    *----*----* |/|
*   |4* | 2  |  3 |    | 4  |  5 |5*3|
*   |/|2*----*----*    |    |    |/| *
*   * |/ 2  /  3 /     *----*----* |/
*   |0*----*----*      |    |    |1*
*   |/0   /  1 /       | 0  |  1 |/
*   *----*----*        *----*----*
* @endverbatim
*
* *---*---*
* | 2 | 3 |
* *---*---*    case_xy      (one isotropic refinement step)
* | 0 | 1 |
* *---*---*
*
* *---*---*
* |   |   |
* | 0 | 1 |    case_x
* |   |   |
* *---*---*
*
* *-------*
* |   1   |
* *-------*    case_y
* |   0   |
* *-------*
*
*/

template <EquationsType equationsType, int dim>
bool AssemblingUtilities<equationsType, dim>::is_refinement_compatible_with_subface(int face_no, RefinementCase<dim> ref_case, int child, int subface_no, RefinementCase<dim - 1> face_ref_case)
{
  if (ref_case == RefinementCase<dim>::cut_xyz)
  {
    switch (face_no) {
    case 0:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 2 && subface_no == 1) return true;
        if (child == 4 && subface_no == 2) return true;
        if (child == 6 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 2 && subface_no == 1) return true;
        if (child == 4 && subface_no == 0) return true;
        if (child == 6 && subface_no == 1) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 2 && subface_no == 0) return true;
        if (child == 4 && subface_no == 1) return true;
        if (child == 6 && subface_no == 1) return true;
      }
      break;
    case 1:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 1 && subface_no == 0) return true;
        if (child == 3 && subface_no == 1) return true;
        if (child == 5 && subface_no == 2) return true;
        if (child == 7 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        if (child == 1 && subface_no == 0) return true;
        if (child == 3 && subface_no == 1) return true;
        if (child == 5 && subface_no == 0) return true;
        if (child == 7 && subface_no == 1) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        if (child == 1 && subface_no == 0) return true;
        if (child == 3 && subface_no == 0) return true;
        if (child == 5 && subface_no == 1) return true;
        if (child == 7 && subface_no == 1) return true;
      }
      break;
    case 2:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 4 && subface_no == 1) return true;
        if (child == 1 && subface_no == 2) return true;
        if (child == 5 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 4 && subface_no == 1) return true;
        if (child == 1 && subface_no == 0) return true;
        if (child == 5 && subface_no == 1) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 4 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
        if (child == 5 && subface_no == 1) return true;
      }
      break;
    case 3:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 2 && subface_no == 0) return true;
        if (child == 6 && subface_no == 1) return true;
        if (child == 3 && subface_no == 2) return true;
        if (child == 7 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        if (child == 2 && subface_no == 0) return true;
        if (child == 6 && subface_no == 1) return true;
        if (child == 3 && subface_no == 0) return true;
        if (child == 7 && subface_no == 1) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        if (child == 2 && subface_no == 0) return true;
        if (child == 6 && subface_no == 0) return true;
        if (child == 3 && subface_no == 1) return true;
        if (child == 7 && subface_no == 1) return true;
      }
      break;
    case 4:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
        if (child == 2 && subface_no == 2) return true;
        if (child == 3 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
        if (child == 2 && subface_no == 0) return true;
        if (child == 3 && subface_no == 1) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 0) return true;
        if (child == 2 && subface_no == 1) return true;
        if (child == 3 && subface_no == 1) return true;
      }
      break;
    case 5:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 4 && subface_no == 0) return true;
        if (child == 5 && subface_no == 1) return true;
        if (child == 6 && subface_no == 2) return true;
        if (child == 7 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        if (child == 4 && subface_no == 0) return true;
        if (child == 5 && subface_no == 1) return true;
        if (child == 6 && subface_no == 0) return true;
        if (child == 7 && subface_no == 1) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        if (child == 4 && subface_no == 0) return true;
        if (child == 5 && subface_no == 0) return true;
        if (child == 6 && subface_no == 1) return true;
        if (child == 7 && subface_no == 1) return true;
      }
      break;
    }
  }
  else if (ref_case == RefinementCase<dim>::cut_xy)
  {
    switch (face_no) {
    case 0:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 2 && subface_no == 1) return true;
        if (child == 0 && subface_no == 2) return true;
        if (child == 2 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 2 && subface_no == 1) return true;
        if (child == 0 && subface_no == 0) return true;
        if (child == 2 && subface_no == 1) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        exit(1);
      }
      break;
    case 1:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 1 && subface_no == 0) return true;
        if (child == 3 && subface_no == 1) return true;
        if (child == 1 && subface_no == 2) return true;
        if (child == 3 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        if (child == 1 && subface_no == 0) return true;
        if (child == 3 && subface_no == 1) return true;
        if (child == 1 && subface_no == 0) return true;
        if (child == 3 && subface_no == 1) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        exit(1);
      }
      break;
    case 2:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 0 && subface_no == 1) return true;
        if (child == 1 && subface_no == 2) return true;
        if (child == 1 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        exit(1);
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
        if (child == 1 && subface_no == 1) return true;
      }
      break;
    case 3:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 2 && subface_no == 0) return true;
        if (child == 2 && subface_no == 1) return true;
        if (child == 3 && subface_no == 2) return true;
        if (child == 3 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        exit(1);
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        if (child == 2 && subface_no == 0) return true;
        if (child == 2 && subface_no == 0) return true;
        if (child == 3 && subface_no == 1) return true;
        if (child == 3 && subface_no == 1) return true;
      }
      break;
    case 4:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
        if (child == 2 && subface_no == 2) return true;
        if (child == 3 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
        if (child == 2 && subface_no == 0) return true;
        if (child == 3 && subface_no == 1) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 0) return true;
        if (child == 2 && subface_no == 1) return true;
        if (child == 3 && subface_no == 1) return true;
      }
      break;
    case 5:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
        if (child == 2 && subface_no == 2) return true;
        if (child == 3 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
        if (child == 2 && subface_no == 0) return true;
        if (child == 3 && subface_no == 1) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 0) return true;
        if (child == 2 && subface_no == 1) return true;
        if (child == 3 && subface_no == 1) return true;
      }
      break;
    }
  }
  else if (ref_case == RefinementCase<dim>::cut_y)
  {
    switch (face_no) {
    case 0:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 0 && subface_no == 2) return true;
        if (child == 1 && subface_no == 1) return true;
        if (child == 1 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        return true;
      }
      break;
    case 1:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 0 && subface_no == 2) return true;
        if (child == 1 && subface_no == 1) return true;
        if (child == 1 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        return true;
      }
      break;
    case 2:
      return true;
      break;
    case 3:
      return true;
      break;
    case 4:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 0 && subface_no == 1) return true;
        if (child == 1 && subface_no == 2) return true;
        if (child == 1 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
      }
      break;
    case 5:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 0 && subface_no == 1) return true;
        if (child == 1 && subface_no == 2) return true;
        if (child == 1 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
      }
      break;
    }
  }

  else if (ref_case == RefinementCase<dim>::cut_x)
  {
    switch (face_no) {
    case 0:
      return true;
      break;
    case 1:
      return true;
      break;
    case 2:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 0 && subface_no == 1) return true;
        if (child == 1 && subface_no == 2) return true;
        if (child == 1 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
      }
      break;
    case 3:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 0 && subface_no == 1) return true;
        if (child == 1 && subface_no == 2) return true;
        if (child == 1 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
      }
      break;
    case 4:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 0 && subface_no == 2) return true;
        if (child == 1 && subface_no == 1) return true;
        if (child == 1 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        return true;
      }
      break;
    case 5:
      if (face_ref_case == RefinementCase<dim - 1>::cut_xy)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 0 && subface_no == 2) return true;
        if (child == 1 && subface_no == 1) return true;
        if (child == 1 && subface_no == 3) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_x)
      {
        if (child == 0 && subface_no == 0) return true;
        if (child == 1 && subface_no == 1) return true;
      }
      else if (face_ref_case == RefinementCase<dim - 1>::cut_y)
      {
        return true;
      }
      break;
    }
  }
}

template <EquationsType equationsType, int dim>
bool AssemblingUtilities<equationsType, dim>::is_refinement_within_current_cell(int face_no, RefinementCase<dim> ref_case, int child)
{
  switch (face_no) {
  case 0:
    if (ref_case == RefinementCase<dim>::cut_xy)
      return ((child == 1) || (child == 3));
    else if (ref_case == RefinementCase<dim>::cut_xyz)
      return ((child == 1) || (child == 3) || (child == 5) || (child == 7));
    else if (ref_case == RefinementCase<dim>::cut_y)
      return false;
    else if (ref_case == RefinementCase<dim>::cut_x)
      return (child == 1);
    else {
      std::cout << "Unknown refinement" << ref_case;
      exit(1);
    }
    break;
  case 1:
    if (ref_case == RefinementCase<dim>::cut_xy)
      return (!((child == 1) || (child == 3)));
    else if (ref_case == RefinementCase<dim>::cut_xyz)
      return (!((child == 1) || (child == 3) || (child == 5) || (child == 7)));
    else if (ref_case == RefinementCase<dim>::cut_y)
      return false;
    else if (ref_case == RefinementCase<dim>::cut_x)
      return (child == 0);
    else {
      std::cout << "Unknown refinement";
      exit(1);
    }
    break;
  case 2:
    if (ref_case == RefinementCase<dim>::cut_xy)
      return ((child == 2) || (child == 3));
    else if (ref_case == RefinementCase<dim>::cut_xyz)
      return ((child == 2) || (child == 3) || (child == 6) || (child == 7));
    else if (ref_case == RefinementCase<dim>::cut_y)
      return (child == 1);
    else if (ref_case == RefinementCase<dim>::cut_x)
      return false;
    else {
      std::cout << "Unknown refinement";
      exit(1);
    }
    break;
  case 3:
    if (ref_case == RefinementCase<dim>::cut_xy)
      return (!((child == 2) || (child == 3)));
    else if (ref_case == RefinementCase<dim>::cut_xyz)
      return (!((child == 2) || (child == 3) || (child == 6) || (child == 7)));
    else if (ref_case == RefinementCase<dim>::cut_y)
      return (child == 0);
    else if (ref_case == RefinementCase<dim>::cut_x)
      return false;
    else {
      std::cout << "Unknown refinement";
      exit(1);
    }
    break;
  case 4:
    if (ref_case == RefinementCase<dim>::cut_xy)
      return false;
    else if (ref_case == RefinementCase<dim>::cut_xyz)
      return (child > 3);
    else if (ref_case == RefinementCase<dim>::cut_y)
      return false;
    else if (ref_case == RefinementCase<dim>::cut_x)
      return false;
    else {
      std::cout << "Unknown refinement";
      exit(1);
    }
    break;
  case 5:
    if (ref_case == RefinementCase<dim>::cut_xy)
      return false;
    else if (ref_case == RefinementCase<dim>::cut_xyz)
      return (child <= 3);
    else if (ref_case == RefinementCase<dim>::cut_y)
      return false;
    else if (ref_case == RefinementCase<dim>::cut_x)
      return false;
    else {
      std::cout << "Unknown refinement";
      exit(1);
    }
    break;
  }
  exit(1);
  return true;
}

template class AssemblingUtilities<EquationsTypeMhd, 3>;
