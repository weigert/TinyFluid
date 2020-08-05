# TinyFluid
Finite-Volume based Navier-Stokes and general Transport-Equation / Conservation-Law solver

Originally written to solve fluid-dynamics problems to model physical-phenomena for procedural world generation.

Written also to be an understandable implementation of the SIMPLE / SIMPLEC algorithm in C++, based on the book "Computational Methods for Fluid Dynamics" by Ferziger and Peric. Also implements the solution of the general transport equations (differential conservation laws / PDEs) for arbitrary sources and sinks, using numerical methods, in a legible manner.

If you have questions or find something unclear, feel free to open an issue.

![Pressure](https://github.com/weigert/TinyFluid/blob/master/screenshots/pressure.png)

*Time-evolution of two gaussian pressure spikes on the central x-axis dissipating and then colliding, forming vortices. Lines indicate direction, color indicates velocity magnitude.*

## Features

- Explicit and implicit time-discretization schemes of higher orders, and corresponding support for point-implicit solutions
- Arbitrary order space-discretization operators using finite differences
- Arbitrary order space-interpolation operators using lagrange polynomials
- Support for 1D, 2D and 3D regular orthogonal grids (solution generated using flattened arrays though)
- Support for arbitrary source terms for coupling transport equations

## Usage
To include the numerical solver in your own program, copy the folder `TinyFluid` to your program and include the main file `TinyFluid.h`. Then you have access to the methods.

See the individual examples (e.g. the base example) for information on how to use the renderer, as well as how to construct and then solve a fluid-mechanics problem.

### Compilation

Compiled using gcc C++14 on Ubuntu 18 LTS

Use the make file included to compile the examples:

    sudo make all

Make sure to properly link the dependencies listed below to compile your own program. See makefile for examples.

Note: Make sure to compile with the `-O3` flag so that Eigen is optimized. This will make it run orders of magnitude faster.

### Dependencies

#### ProcFluid Core:
    - Eigen C++ (Linear Algebra)
    - GLM (Useful small-vector math)

#### Renderer:
    - SDL2 (Optional)

#### Other
    - Timer Helper: pthread

## To-Do

There are still some issues with the pressure correction laplace equation concerning spurious modes appearing and causing checkerboard patterns. I need to read up again on why exactly there are problems with choosing an FFD + BFD scheme. It's hard to wrap my head around.


    - Arbitrary order space FD
    - Arbitrary order space interpolation
    - Simple transporter function for a known flow-field
    - Write a Wiki that contains the necessary information on how this works
        - Time Integration
        - Space Discretizations
        - Boundary and Initial Conditions
        - Example Problem Closer Explanation
        - Included Source Terms
        - Renderer Breakdown
        - Helpers
    - Make the solver capable of handling variable cell volumes and contact areas
    - Make the solver *theoretically* capable of handling 3D systems of cells

    - Write an example for ocean-current computation
    - Write an example for ground-water flow computation

## Sources
This project was based on information from the book "Computational Methods for Fluid Dynamics" by Ferziger and Peric. Particularly chapter 7 of the latest edition (at the time of writing - January 2020) helped me in implementing the solution of the pressure coupled Navier-Stokes equation.

The code is still written 100% from scratch, as I found the Fortran77 implementation referenced in the book absolutely unintelligible.

A big thank you to Dr. Daniel Meyer-Massetti for his excellent course "Numerical Methods for Transport Phenomena" at ETH Zurich that first introduced me to numerical methods for solving PDEs. The content of that course was the strong foundation that even let me know how to begin to create this repository.

## License
All source code is licensed under the MIT License.
