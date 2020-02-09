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
  void solveNavierStokes(Eigen::VectorXd& vX, Eigen::VectorXd& vY, Eigen::VectorXd& P, Eigen::SparseMatrix<double> B, Eigen::SparseMatrix<double> GXF,  Eigen::SparseMatrix<double> GXB,  Eigen::SparseMatrix<double> GYF,  Eigen::SparseMatrix<double> GYB, double tol, int maxiter, double relaxP, double relaxV){
    std::cout<<"Please specifiy a valid solution method."<<std::endl;
  }

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;

  template<> //Naive Solution (Strange Pressure Behavior)
  void solveNavierStokes<SIMPLEC>(Eigen::VectorXd& vX, Eigen::VectorXd& vY, Eigen::VectorXd& P, Eigen::SparseMatrix<double> B, Eigen::SparseMatrix<double> GXF,  Eigen::SparseMatrix<double> GXB,  Eigen::SparseMatrix<double> GYF,  Eigen::SparseMatrix<double> GYB, double tol, int maxiter, double relaxP, double relaxV){
    solver.compute(GXF*GXB + GYF*GYB);
    solver.setTolerance(0.000001);


    //
    Eigen::VectorXd dP;
    Eigen::VectorXd dvX;
    Eigen::VectorXd dvY;

    double newerr = 1.0;
    bool divergence = false;
    while(newerr > tol && !divergence && maxiter){
      //Pressure Correction

      dP = solver.solve(GXB*vX + GYB*vY);

      // We compute the velocity correction based on the pressure correction

      dvX = -(GXF*dP);
      dvY = -(GYF*dP);

      // We correct our velocity guesses from the intermediary field

      vX += relaxV*B*dvX;
      vY += relaxV*B*dvY;
      P += relaxP*B*dP;

      // Compute the Error (Divergence of Field) simply using our guess.
      newerr = (GXB*vX + GYB*vY).squaredNorm()/(SIZE*SIZE);
      maxiter--;

      //std::cout<<newerr<<std::endl;
      if(newerr > 1E+5) divergence = true;
    }

    //Error Handling
    if(!maxiter) std::cout<<"Max iterations surpassed."<<std::endl;
    if(divergence) std::cout<<"Instability encountered."<<std::endl;
  };

  /*
    Transport Equation Integrator:

      Define your own parameter set...

  */

  void solveTransport(double dt, double dx, double dy, Eigen::VectorXd& val, Eigen::VectorXd& vX, Eigen::VectorXd& vY){
    double diffusivity = 0.001;

    //Velocity Matrices
    Eigen::SparseMatrix<double> vX_MAT = alg::sparseDiagonalize(vX);
    Eigen::SparseMatrix<double> vY_MAT = alg::sparseDiagonalize(vY);

    Eigen::SparseMatrix<double> XFLUX = space::FV_FLUX(glm::vec2(1, 0))/dx;
    Eigen::SparseMatrix<double> YFLUX = space::FV_FLUX(glm::vec2(0, 1))/dy;

    Eigen::SparseMatrix<double> XDIFFUSION = space::FV_DIFFUSION(glm::vec2(1, 0))/dx/dx;
    Eigen::SparseMatrix<double> YDIFFUSION = space::FV_DIFFUSION(glm::vec2(0, 1))/dy/dy;

    //Get our Operators
    Eigen::SparseMatrix<double> CONVECTION = -vX_MAT*XFLUX - vY_MAT*YFLUX;
    Eigen::SparseMatrix<double> DIFFUSION = diffusivity*(XDIFFUSION + YDIFFUSION);

    //Full Matrix Operator...
    Eigen::SparseMatrix<double> MAT = CONVECTION + DIFFUSION;

    //Integrate!
    PDE::integrate<PDE::CN>(dt, val, MAT);
  }
};
