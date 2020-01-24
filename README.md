# PDE-Transport
Finite-Volume based Navier-Stokes and General Transport Equation / Conservation Law Solver

This repository contains the code for solving abritrary transport equations (i.e. differential conservation laws), or PDEs with those terms, using numerical methods.

It provides a simple interface with a number of time-integrators to choose from, as well as different order approximations for various terms.

You can also choose to include source terms to simply couple your equations.

The system will also solve the Navier-Stokes equations if the underlying flow-field is unknown. It does this using a pressure correction method.

This uses Eigen in C++ for the linear algebra.

This repository has more of a teaching purpose, so if you are interested in how to do your own CFD simulations, the code is written to be legible and comprehensible. If you find something unclear, please contact me.

## License
MIT License
