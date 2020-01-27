/*

PDE Transport Master File

  Currently, independent of the time-step, we are accumulating error at the corners for some reason.
  This needs to be fixed.



My operator discretizations are fine...

You could also do weighted sums of different types of approximations! (e.g. upwind + central finite differences)
I should implement an upwind operator too.

The pressure field we want to use for the new velocity field should be the one at the next step.
We don't know this one though.

We compute a new velocity field from an imperfect pressure field (i.e. not the correct one).

We must compute the pressure correction value
Then we must compute the speed correction value


1. Have some initial guess P* for pressure field (previous field)
2. Solve Impulse Equations for U*, V*
3. Solve the pressure correction equation for dP, get dU and dV
4. Correct P = P* + dP, U = U* + dU, V = V* + dV
5. Set P* = P, Repeat from Step 2

Pressure correction can converge if we don't underrelax.

*/

//Eigen Stuff
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <glm/glm.hpp>

#include <iostream>

// Core Stuff
#include "src/spacesolve.h"   //Space Discretizations
#include "src/timesolve.h"    //Time Integrators
#include "src/fullsolve.h"    //Full Transport / NS Solvers

// Rendering Stuff
#include "render/view.cpp"    //Renderer (Requires SDL)
#include "render/input.cpp"    //Renderer (Requires SDL)
