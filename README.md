# PDE-Transport
Finite-Volume based Navier-Stokes and General Transport Equation / Conservation Law Solver

This repository contains the code for solving abritrary transport equations (i.e. differential conservation laws), or PDEs with those terms, using numerical methods.

It provides a simple interface with a number of time-integrators to choose from, as well as different order approximations for various terms.

You can also choose to include source terms to simply couple your equations.

The system will also solve the Navier-Stokes equations if the underlying flow-field is unknown. It does this using a pressure correction method.

This uses Eigen in C++ for the linear algebra.

This repository has more of a teaching purpose, so if you are interested in how to do your own CFD simulations, the code is written to be legible and comprehensible. If you find something unclear, please contact me.

## Applications
This allows you to solve a number of problems and some examples have been included.

The initial purpose is to solve the Navier Stokes equations on a topographical height-map to get the wind-patterns, and to use the resulting velocity field to solve the transport equations for humidity and temperature.

By introducing source terms that couple these (e.g. Raoult's Law for the phase of the water, Henry's Law to derive mass-transfer into the air over bodies of water, etc.) we can extract dynamic weather effects.

## Usage
This has been written in library form to make it easy to include the PDE solver in your application.

A rendering helper has been included for practicality, if you wish to write a test program and visualize the information quickly.

### Dependencies

- Eigen C++

## To-Do
- Higher orders (automatically) for the various discrete operators
- Simple transporter function for a known flow-field
- Write a Wiki that contains the necessary information on how this works
- Make the solver capable of handling variable cell volumes and contact areas
- Make the solver *theoretically* capable of handling 3D systems of cells
- 

- Write an example for ocean-current computation
- Write an example for ground-water flow computation

## License
MIT License
