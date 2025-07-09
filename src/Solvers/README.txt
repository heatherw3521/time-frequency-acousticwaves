This folder consists of the following functions:

`solve_HH_density`
The density layer is iteratively, for each frequency k, solved through `HH_density`.

`HH_density`
The density layer is computed for a single frequency through Zetatrap. Other methods to solve for the density layer given one frequency will be here.

`solve_HH_eqn`
Evaluates the solution to the Helmholtz equation using the density layer gotten from `solve_HH_desnity` through `eval_sol_freqs` to be used in `solve_HH_to_wave`.

`eval_sol_freqs`
Evaluates the solution to the Helmholtz equation on given points iteratively, for each frequency k, through zetatap or fmm.

`solve_HH_to_wave`
Evaluates the solution to the wave equation using information gotten from `solve_HH_eqn` using the complex deformation, fast transform or Gauss-legendre quadrature through `eval_sol_time`. (wrapper)

`eval_sol_time`
Evaluates the solution to the wave equation at domain points using using the complex deformation, fast transform or Gauss-legendre quadrature.