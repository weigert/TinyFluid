/*
  Full Navier-Stokes / Transport Equation Solvers!
*/

namespace PDE{

  enum Solver{
    TRANSPORT,
    NAIVE,
    CHOI_NONIT
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
  void solveNavierStokes<NAIVE>(double dt, Eigen::VectorXd& vX, Eigen::VectorXd& vY){
    /*
    //Initialize Volume and Area Maps
    P = P0*Eigen::ArrayXd::Ones(SIZE*SIZE);

    //Initialize Velocities to Zero
    vX = u*Eigen::ArrayXd::Ones(SIZE*SIZE);
    vY = v*Eigen::ArrayXd::Ones(SIZE*SIZE);
    vX += u*am::flatGaussian();
    vY += v*am::flatGaussian();

    vX_MAT = am::sparseDiagonalize(vX);
    vY_MAT = am::sparseDiagonalize(vY);

    //Get the Transport Operator
    fX_f = am::getTrapezeX(1);
    fX_b = am::getTrapezeX(-1);
    fY_f = am::getTrapezeY(1);
    fY_b = am::getTrapezeY(-1);

    dX_f = am::getDiffX(1);
    dX_b = am::getDiffX(-1);
    dY_f = am::getDiffY(1);
    dY_b = am::getDiffY(-1);

    //Discrete Differential Operators (of some order)
    LX = am::laplaceX();
    LY = am::laplaceY();
    GX = am::gradX();
    GY = am::gradY();


    //Loop Here


    Eigen::SparseMatrix<double> TXX = viscosity * (dX_f/dx/dx + dX_b/dx/dx);
    Eigen::SparseMatrix<double> TXY = viscosity * (dY_f/dy/dy + dY_b/dy/dy);
    Eigen::SparseMatrix<double> TYX = viscosity * (dX_f/dx/dx + dX_b/dx/dx);
    Eigen::SparseMatrix<double> TYY = viscosity * (dY_f/dy/dy + dY_b/dy/dy);

    vX_MAT = am::sparseDiagonalize(vX);
    vY_MAT = am::sparseDiagonalize(vY);

    J = -vY_MAT*((fY_f-fY_b)/dy)-vX_MAT*((fX_f-fX_b)/dx);   //This includes the velocity field
    
    Eigen::SparseMatrix<double> MATX = J+TXX+TXY;
    Eigen::SparseMatrix<double> MATY = J+TYX+TYY;


    double dt = 0.0001;

    Eigen::VectorXd vX_Temp = vX;
    Eigen::VectorXd vY_Temp = vY;

    //First Explicit Integration Step (Without Pressure Source)
    PDE::integrate<PDE::CN>(dt, vX_Temp, MATX);
    PDE::integrate<PDE::CN>(dt, vY_Temp, MATY);

    //Compute Pressure Field EXPLICITLY (Laplace Equation)
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute((LX/dx/dx + LY/dy/dy));                 //2D Laplace Operator
    P = solver.solve((GX/dx*vX_Temp + GY/dy*vY_Temp)/dt);  //Maybe this dt needs to go
    std::cout<<P<<std::endl;

    //Explicit
    vX = vX_Temp - dt*((fX_f-fX_b)/dx*P);
    vY = vY_Temp - dt*((fY_f-fY_b)/dy*P);
    */
  };


  template<> //Choi et. al (1994) (NON-ITERATIVE, req. small dt)
  void solveNavierStokes<CHOI_NONIT>(double dt, Eigen::VectorXd& vX, Eigen::VectorXd& vY){

    /*
    Missing: The Operators!

    MATX
    MATY
    dx
    dy

    */
/*
    //Compute the Pressure Source-Terms for the Initial Guess
    Eigen::VectorXd P_SOURCE_X = -2.0*((fX_f-fX_b)/dx)*P;
    Eigen::VectorXd P_SOURCE_Y = -2.0*((fY_f-fY_b)/dy)*P;

    //Storage for Pressure Correction
    Eigen::VectorXd dP;

    //Temporary
    Eigen::VectorXd vX_Temp = vX;
    Eigen::VectorXd vY_Temp = vY;

    //Sparse Eigen-Solver
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute((LX/dx/dx + LY/dy/dy));                      //2D Laplace Operator

    //This section could be converged! (To get a better in-between estimate)
      //Crank-Nicholson Integration (With full Explicit Source)
      PDE::integrate<PDE::EE>(0.5*dt, vX_Temp, MATX, P_SOURCE_X);
      PDE::integrate<PDE::IE>(0.5*dt, vX_Temp, MATX);
      PDE::integrate<PDE::EE>(0.5*dt, vY_Temp, MATY, P_SOURCE_Y);
      PDE::integrate<PDE::IE>(0.5*dt, vY_Temp, MATY);

      //Pressure Correction
      dP = solver.solve((GX/dx*vX_Temp + GY/dy*vY_Temp)*2.0/dt);       //Maybe this dt needs to go

    //Final Calculation
    std::cout<<"Final Calculation..."<<std::endl;
    P_SOURCE_X = -((fX_f-fX_b)/dx)*dP;        //New Values from Pressure Correction
    P_SOURCE_Y = -((fY_f-fY_b)/dy)*dP;
    vX = vX_Temp + dt * 0.5 * P_SOURCE_X;     //Only Consider Pressure Correction!!
    vY = vY_Temp + dt * 0.5 * P_SOURCE_Y;

    //Fix the Pressure
    P += dP;
    */
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
