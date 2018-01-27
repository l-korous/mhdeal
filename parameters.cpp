#include "parameters.h"
#include "equationsMhd.h"
#include "DealiiExtensions.h"

Parameters<>::Parameters()
{
  parameters.initCond = 0;
  parameters.num_flux_type = hlld;
  parameters.initial_and_max_cfl_constant = parameters.cfl_constant = .05;
  parameters.quadrature_order = 5;
  parameters.initial_quadrature_order = 10;
  parameters.polynomial_order_dg = 1;
  parameters.polynomial_order_hdiv = 0;
  parameters.limit_in_nonlin_loop = true;
  
  
  parameters.initial_and_max_newton_damping = .9;
  parameters.decrease_factor = .8;
  parameters.increase_factor = 1.1;

  parameters.patches = 2;
  parameters.output_step = -1.e-3;

  parameters.debug = false;
  parameters.debug_limiter = false;
  parameters.debug_dofs = false;

  parameters.output_matrix = false;
  parameters.output = OutputType::quiet_solver;
  parameters.output_rhs = false;
  parameters.output_solution = false;

  parameters.snapshot_step = 1.;

  parameters.time_step = 1.e-5;
  parameters.final_time = 10.;

  parameters.solver = gmres;
  parameters.linear_residual = 1e-10;
  parameters.max_iterations = 10000;
  parameters.ilut_fill = 1.5;
  parameters.ilut_drop = 1e-6;
  parameters.ilut_atol = 1e-6;
  parameters.ilut_rtol = 1.0;

  parameters.gas_gamma = 1.4;

  parameters.newton_max_iterations = 100;
  parameters.newton_residual_norm_threshold = 1e-8;

  parameters.lax_friedrich_stabilization_value = 0.;

  parameters.is_stationary = false;

  parameters.needs_gradients = false;
  parameters.needs_forcing = false;
}

template class Parameters<2>;
template class Parameters<3>;