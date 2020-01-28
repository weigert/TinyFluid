/*
  Full Navier-Stokes / Transport Equation Solvers!
*/

namespace PDE{

  enum Solver{
    TRANSPORT,
    SIMPLEC
  };

  //Full-NS-Solver
  template<Solver S>
  void solveNavierStokes(double dt, Eigen::VectorXd& vX, Eigen::VectorXd& vY){
    std::cout<<"Please specifiy a valid solution method."<<std::endl;
  }

  //Full Transport Solver
  template<Solver S>
  void solveTransport(double dt, Eigen::VectorXd& val, Eigen::VectorXd& vX, Eigen::VectorXd& vY, Eigen::VectorXd& source){
    std::cout<<"Please specifiy a valid solution method."<<std::endl;
  }

  template<> //Naive Solution (Strange Pressure Behavior)
  void solveNavierStokes<SIMPLEC>(double dt, Eigen::VectorXd& vX, Eigen::VectorXd& vY){
    return;
  };

  /* Simple Transport Integrator! */
  template<>
  void solveTransport<TRANSPORT>(double dt, Eigen::VectorXd& val, Eigen::VectorXd& vX, Eigen::VectorXd& vY, Eigen::VectorXd& source){

    /*
    Get the appropriate transport operators...
    Then include the source term.

    We include the velocity field!

    And the area is important too actually.

    //DIFFUSION OPERATOR
    G = viscosity*((dY_f+dY_b)/dy/dy + (dX_f+dX_b)/dx/dx);

    //Convection Operator!
    J = -vY_MAT*((fY_f-fY_b)/dy)-vX_MAT*((fX_f-fX_b)/dx);

    //Full Operator
    MAT = J+G;

    //PDE::integrate<PDE::EE>(0.0005, heightMap, MAT);
    //PDE::integrate<PDE::IE>(0.01, heightMap, MAT);
    PDE::integrate<PDE::CN>(0.01, heightMap, MAT);
    */

    return;
  }
};
