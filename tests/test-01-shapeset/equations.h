#ifndef _EQUATIONS_H
#define _EQUATIONS_H

enum EquationsType {
  EquationsTypeEuler,
  EquationsTypeMhd
};

// Dummy class for templating.
template <int dim>
class Equations
{
public:
  // Self-explanatory
  static const unsigned int n_components = 2 * dim + 2;
  static const unsigned int first_momentum_component = 1;
  static const unsigned int first_magnetic_flux_component = dim + 1;
  static const unsigned int density_component = 0;
  static const unsigned int energy_component = 2 * dim + 1;

  static std::vector<std::string> component_names();
  static std::vector<DataComponentInterpretation::DataComponentInterpretation> component_interpretation();
};

#endif