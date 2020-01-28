# ProcFluid
Finite-Volume based Navier-Stokes and general Transport-Equation / Conservation-Law solver

Originally written to solve fluid-dynamics problems to model physical-phenomena for procedural world generation.

Blog article for this repository and why I made it: [Fluid-Dynamics for Procedural Climate Systems](weigert.vsos.ethz.ch)

The **Wiki** for this project contains more information about the features, implementation, where to find more information on the various subjects, as well as the explanations of a number of example problems solved using this library.

Written also to be an understandable implementation of the SIMPLE / SIMPLEC algorithm in C++, based on the book "Computational Methods for Fluid Dynamics" by Ferziger and Peric. Also implements the solution of the general transport equations (differential conservation laws / PDEs) for arbitrary sources and sinks, using numerical methods, in a legible manner.

If you have questions or find something unclear, feel free to open an issue.

![Pressure](https://github.com/weigert/ProcFluid/blob/master/screenshots/pressure.png)

*Time-evolution of two gaussian pressure spikes on the central x-axis dissipating and then colliding, forming vortices. Lines indicate direction, color indicates velocity magnitude.*


## Features
### CFD / Numerical Solver
    - Simple, understandable interface for defining and solving a fluid-mechanics problem
    - Fully visible implementation of all algorithms for adaptation / implementation of other solvers
    - Explicit and implicit time-discretization schemes of higher orders, and corresponding support for point-implicit solutions
    - Arbitrary order space-discretization operators using finite differences
    - Arbitrary order space-interpolation operators using lagrange polynomials
    - Support for 1D, 2D and 3D regular orthogonal grids (solution generated using flattened arrays though)
    - Support for arbitrary source terms for coupling equations; some physical equations are given as examples

### Visualizer
    - Simple rendering helper class that creates a window and then interfaces like an HTML canvas for rendering
    - Very bare-bones, easy to understand implementation with SDL2
    - Drawing helper methods let you visualize different data in different ways
    - Arbitrary color-schemes using bezier curve interpolation

Feel free to steal the renderer for your C++ application if you need quick visualization of 2D information without writing a whole visualization wrapper. Can be included entirely independently.

### Other
    - Timer helper class for benchmarking arbitrary code / executing code at fixed intervals in a detached thread

## Examples / Applications
This library allows you to solve a number of problems and some examples have been included below.

### Pressure Wave and Velocity Field (Navier-Stokes Base Example)
... describe here and add some images...

*Folder:* `examples/flowtest`

### Topography and Climate Simulation
We can solve the pressure-linked Navier-Stokes equations on a topographical height-ma in order to generate wind-patterns. We can do so by approximating it as a 2D problem, where the volume of individual cells is proportional to (1-height), if the heightmap is scaled between \[0, 1\]. This yields the pressure and velocity fields (i.e. wind-patterns).

Subsequently, other quantities (temperature, humidity) can simple be transported by the velocity field with defined diffusivities.

Coupling the two equations using physical source terms (i.e. raoults law for water phase, henrys law for mass-transfer from bodies of water, and a defined heat of evaporation), we can simulate dynamic weather effects.

### Ocean Current Simulation

## Usage
To include the numerical solver in your own program, copy the folder pde-transport to your program and include the main file `pde-transport.h`. Then you have access to the methods.

See the individual examples (e.g. the base example) for information on how to use the renderer, as well as how to construct and then solve a fluid-mechanics problem.

### Compilation

Compiled using gcc C++14 on Ubuntu 18 LTS

Use the make file included to compile the examples:

    sudo make all
    
Make sure to properly link the dependencies listed below to compile your own program. See makefile for examples.

Note: Make sure to compile with the `-O3` flag so that Eigen is optimized. This will make it run orders of magnitude faster.

### Dependencies

#### PDE-Transport Core:
    - Eigen C++ (Linear Algebra)
    - GLM (Useful small-vector math) 

#### Renderer:
    - SDL2

#### Other
    - Timer Helper: pthread

## To-Do
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
