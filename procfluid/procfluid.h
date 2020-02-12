/*

PDE Transport Master File

*/

//Eigen Stuff
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <glm/glm.hpp>

#include <iostream>

// Core Stuff
#include "src/algebra.h"    //Space Discretizations
#include "src/shape.h"      //Shape Functions

#include "src/space.h"      //Space Discretizations
#include "src/time.h"       //Time Discretizations and Integrators
#include "src/solve.h"      //Full Scheme Solvers / Algorithms
#include "src/source.h"     //Physical Source and Sink Terms

// Rendering Stuff
#include "render/view.cpp"    //Renderer (Requires SDL)
#include "render/input.cpp"    //Renderer (Requires SDL)

#include "helpers/timer.h"    //Benchmarking Tool
