This document describes all that is necessary to insert a new problem so that it can get computed. For the time being, it must be a full MHD problem (in future, there might be an option to also get an Euler problem solved - currently the support is not complete).

1. Equations
  * equations in https://github.com/l-korous/mhdeal/blob/master/equationsMhd.h should be changed when
    * either additional term is needed in the flux, then the method Equations<EquationsTypeMhd, dim>::compute_flux_matrix is what needs to change
   * or additional term that depends on state variables gradients is needed, then the method Equations<EquationsTypeMhd, dim>::compute_jacobian_addition is what needs to change
      * also, in parameters.cpp, you need to set this->needs_gradients = true;
   * or just a source term (term that depends on state variables values) is needed, then Equations<EquationsTypeMhd, dim>::compute_forcing_vector is your friend
      * also, in parameters.cpp, you need to set this->needs_forcing = true;
   * lastly, if you wish to implement a new numerical flux, that would be in Equations<EquationsTypeMhd, dim>::numerical_normal_flux, but you also would need to:
      * add a new entry in Parameters::NumFluxType enum
      * use the new entry in parameters.cpp: this->num_flux_type = <...>
   * there are also relevant parameters like gas_gamma, is_stationary (do you want to solve a stationary or transient problem), theta (as in the theta schema - which is an explicit Euler for theta = 0, Crank-Nicolson for theta = 0.5 and implicit Euler for theta = 1)

2. Triangulation & computation parameters
  * completely delegated to the Parameters class (parameters.h, parameters.cpp)

3. Initial conditions
  * see initialCondition.cpp, namely InitialCondition<EquationsTypeMhd, 3>::value

4. Boundary conditions
  * see boundaryConditions.cpp, namely BoundaryConditions<equationsType, dim>::bc_vector_value
  * The enumeration BoundaryKind (which is defined for a particular set of equaions, e.g. now in https://github.com/l-korous/mhdeal/blob/master/equationsMhd.h#L22) should be used in BoundaryConditions<EquationsTypeMhd, 3>::BoundaryConditions() constructor to assign a particular boundary part the appropriate type.
   * How to obtain boundary_ids from the created triangulation - consult deal.II documentation
   * As we do DG, there is no such thing as Dirichlet conditions, only 'inlet' / 'input' / 'inflow' conditions, whose values get into the solution via the boundary - using numerical fluxes. Therefore the condition that corresponds to Dirichlet is called 'inflow'.
