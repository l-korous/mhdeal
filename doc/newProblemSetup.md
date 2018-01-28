This document describes all that is necessary to insert a new problem so that it can get computed. For the time being, it must be a full MHD problem (in future, there might be an option to also get an Euler problem solved - currently the support is not complete).

0. Running existing examples
  * CMake configuration by default prepares the library (project 'mhdeal') together with several examples (in the 'examples' directory).
  * Build either the entire set of projects, or just the library and the example you are interested in.

0*. Creating your own example
  * The best way is definitely to start off from an existing example, and just fine tune the parameters, triangulation, and possibly equations and / or initial / boundary conditions (see further).

1. Equations
  * equations in https://github.com/l-korous/mhdeal/blob/master/equationsMhd.h should be changed when
    * either additional term is needed in the flux, then the method Equations<EquationsTypeMhd, dim>::compute_flux_matrix is what needs to change
   * lastly, if you wish to implement a new numerical flux, that would be in Equations<EquationsTypeMhd, dim>::numerical_normal_flux, but you also would need to:
      * add a new entry in Parameters::NumFluxType enum
      * use the new entry in parameters.cpp: this->num_flux_type = <...>
   * there are also relevant parameters like gas_gamma, is_stationary (do you want to solve a stationary or transient problem), theta (as in the theta schema - which is an explicit Euler for theta = 0, Crank-Nicolson for theta = 0.5 and implicit Euler for theta = 1)
  * equations in https://github.com/l-korous/mhdeal/blob/master/problem.cpp should be changed when
    * any changes to the equations that are not inside of either the MHD flux, or the numerical flux are needed.
  
2. Triangulation & computation parameters
  * completely delegated to the Parameters class (parameters.h - parameter list, in your main.cpp - parameter values)

3. Initial conditions
  * see initialCondition.cpp, namely InitialCondition<EquationsTypeMhd, 3>::value

4. Boundary conditions
  * see boundaryConditions.cpp, namely BoundaryConditions<equationsType, dim>::bc_vector_value
  * The enumeration BoundaryKind (which is defined for a particular set of equaions, e.g. now in https://github.com/l-korous/mhdeal/blob/master/equationsMhd.h#L22) should be used in BoundaryConditions<EquationsTypeMhd, 3>::BoundaryConditions() constructor to assign a particular boundary part the appropriate type.
   * How to obtain boundary_ids from the created triangulation - consult deal.II documentation
   * As we do DG, there is no such thing as Dirichlet conditions, only 'inlet' / 'input' / 'inflow' conditions, whose values get into the solution via the boundary - using numerical fluxes. Therefore the condition that corresponds to Dirichlet is called 'inflow'.
