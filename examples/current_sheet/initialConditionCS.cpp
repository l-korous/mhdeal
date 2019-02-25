#include "equationsMhd.h"
#include "completeEllipticIntegrals.h"
#include "initialConditionCS.h"

template <EquationsType equationsType, int dim>
InitialConditionCS<equationsType, dim>::InitialConditionCS(Parameters<dim>& parameters, CSParameters& cs_pars) :
      InitialCondition<equationsType, dim>(parameters), cs_parameters(cs_pars)
{
}

template <EquationsType equationsType, int dim>
void InitialConditionCS<equationsType, dim>::vector_value(const std::vector<Point<dim> > &points, std::vector<std::array<double, Equations<equationsType, dim>::n_components> >& value_list) const
{

	for (unsigned int i = 0; i < points.size(); ++i)
	{
		const Point<dim>& p = points[i];
		value_list[i][0] = 1.; //density

		value_list[i][1] = 0.; //momentum density
		value_list[i][2] = 0.;
		value_list[i][3] = 0.;

		value_list[i][5] = 0.0; //magnetic field
		value_list[i][6] = std::tanh(p[0]);
		value_list[i][7] = 0.;
		
	}
	for (unsigned int i = 0; i < points.size(); ++i)
	{
		const Point<dim>& p = points[i];
		value_list[i][4] = (1.15+cs_parameters.beta-Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(value_list[i]))/(2./3) + Equations<EquationsTypeMhd, dim >::compute_magnetic_energy(value_list[i]); //energy density
	}
}

template class InitialConditionCS<EquationsTypeMhd, 3>;
