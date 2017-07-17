#include "util.h"
#include "equations.h"

template <int dim>
std::vector<std::string> Equations<dim>::component_names()
{
  return{ "solution" };
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Equations<dim>::component_interpretation()
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation;
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  return data_component_interpretation;
}

template class Equations<3>;
