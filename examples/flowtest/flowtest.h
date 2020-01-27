#include <noise/noise.h>

class Field{
public:
  //Grid-Size
  double dx = 1.0/(double)SIZE;
  double dy = 1.0/(double)SIZE;
  double dt = 0.001;

  //Parameters
  double viscosity = 0.000148;  //[m^2/s]  //Future: Replace with Temperature Dependent?
  double density = 1.225;       //[kg/m^3]Future: Replace with Temperature Dependent?
  double P0 = 1E+5;             //Initial Pressure [Pa]

  //Volume Force
  double gx = 0.0;
  double gy = 0.0;

  //Integration Stuff
  void initialize();
  void timestep();
  int count = 0;
  bool divergence = false;
  double err = 100.0; //Initial Mass Error

  //Velocities (Actually)
  Eigen::VectorXd vX;
  Eigen::VectorXd vY;
  Eigen::VectorXd P;

  //Laplace Operator
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  Eigen::SparseMatrix<double> LAPLACE;  //Full Laplace Operator in 2D

  //Gradient
  Eigen::SparseMatrix<double> GX;       //Gradient Operator X-Element
  Eigen::SparseMatrix<double> GY;       //Gradient Operator Y-Element

  //Spatial Discretization Operator
  Eigen::SparseMatrix<double> XFLUX;         //Integrated Surface Flux X, UNIFORM AREA!!
  Eigen::SparseMatrix<double> YFLUX;
  Eigen::SparseMatrix<double> XDIFFUSION;    //Integrated Surface Flux X, UNIFORM AREA!!
  Eigen::SparseMatrix<double> YDIFFUSION;
};

void Field::initialize(){
  //Initialize Pressure to some pressures pike in the center...
  P = P0*Eigen::ArrayXd::Ones(SIZE*SIZE);
  P += 100*am::flatGaussian(glm::vec2(SIZE/2.0, SIZE/2.0), SIZE/2.0);
  //P += 100*am::flatGaussian(glm::vec2(2.0*SIZE/3.0, SIZE/2.0), SIZE/3.0);

  //Initialize Velocities to Zero
  vX = Eigen::ArrayXd::Zero(SIZE*SIZE);
  vY = Eigen::ArrayXd::Zero(SIZE*SIZE);

  //Discrete Differential Operator (Laplace)
  LAPLACE = am::laplaceX()/dx/dx + am::laplaceY()/dy/dy;
  solver.analyzePattern(LAPLACE);
  solver.setTolerance(0.000001);

  //Discrete Differential Operators (of some order)
  GX = am::gradX()/dx;  //Gradient needs to be divided by the grid-spacing
  GY = am::gradY()/dy;  //-``-

  //Surface Integrators (i.e. convective fluxes over surface)
  XFLUX = (am::getTrapezeX(1) - am::getTrapezeX(-1))/dx;
  YFLUX = (am::getTrapezeY(1) - am::getTrapezeY(-1))/dy;

  //Diffusion Operator! (Note: This is in FV Form and Includes the Area / Volume!!!)
  XDIFFUSION = (am::getDiffX(1)+am::getDiffX(-1))/dx/dx;
  YDIFFUSION = (am::getDiffY(1)+am::getDiffY(-1))/dy/dy;
}

void Field::timestep(){

  /* SIMPLEC ALGORITHM */

  //Construct Non-Linear Operator from Velocity Field and Other Contributions
  Eigen::SparseMatrix<double> vX_MAT = am::sparseDiagonalize(vX);
  Eigen::SparseMatrix<double> vY_MAT = am::sparseDiagonalize(vY);

  //Implicit and Explicitly Evaluated Operators... (Here: Crank-Nicholson)
  Eigen::SparseMatrix<double> E_MAT_X = -vX_MAT*XFLUX + viscosity*XDIFFUSION;
  Eigen::SparseMatrix<double> E_MAT_Y = -vY_MAT*YFLUX + viscosity*YDIFFUSION;
  Eigen::SparseMatrix<double> E_MAT = 1.0*(E_MAT_X+E_MAT_Y);

  Eigen::SparseMatrix<double> I_MAT_X = 0.0*E_MAT_X;
  Eigen::SparseMatrix<double> I_MAT_Y = 0.0*E_MAT_Y;
  Eigen::SparseMatrix<double> I_MAT = 0.0*E_MAT;

  //Construct the weight vector??? ONLY FROM THE IMPLICIT MATRIX PORTION!!
  Eigen::VectorXd E = Eigen::ArrayXd::Ones(SIZE*SIZE);

  //All Elements! (Summed for Each Individual Guy)
  Eigen::VectorXd A_X = ((am::sparseIdentity()-dt*I_MAT_X)*E);
  Eigen::VectorXd A_Y = ((am::sparseIdentity()-dt*I_MAT_Y)*E);

  //All Elements!
  Eigen::VectorXd UA_X = E.cwiseQuotient(A_X);
  Eigen::VectorXd UA_Y = E.cwiseQuotient(A_Y);

  //Refactorize the Solver Matrix
  //solver.factorize(am::sparseDiagonalize(UA_X)*am::laplaceX()/dx/dx + am::sparseDiagonalize(UA_Y)*am::laplaceY()/dy/dy); //2D Laplace Operator
  solver.factorize(GX*am::sparseDiagonalize(UA_X)*GX + GY*am::sparseDiagonalize(UA_Y)*GY); //2D Laplace Operator

  Eigen::VectorXd P_SOURCE_X;
  Eigen::VectorXd P_SOURCE_Y;

  Eigen::VectorXd dP;
  Eigen::VectorXd dvX;
  Eigen::VectorXd dvY;

  //Convergence Iteration
  for(int i = 0; i < 10 && !divergence; i++){

    //Pressure Source-Term from Previous Iteration (and volume forces)
    P_SOURCE_X = gx*Eigen::ArrayXd::Ones(SIZE*SIZE);
    P_SOURCE_Y = gy*Eigen::ArrayXd::Ones(SIZE*SIZE);
    P_SOURCE_X.noalias() += -XFLUX*P;
    P_SOURCE_Y.noalias() += -YFLUX*P;

    //General Integration Scheme
    PDE::integrate<PDE::EE>(dt, vX, E_MAT, P_SOURCE_X);
    PDE::integrate<PDE::EE>(dt, vY, E_MAT, P_SOURCE_Y);
    //PDE::integrate<PDE::IE>(dt, vX, I_MAT);  // MAT could theoretically be split into an implicit and explicit portion
    //PDE::integrate<PDE::IE>(dt, vY, I_MAT);

    //Pressure and Velocity Correction (Intermediate)
    dP.noalias() = solver.solve(GX*vX + GY*vY);
    dvX = -(GX*dP).cwiseProduct(UA_X);
    dvY = -(GY*dP).cwiseProduct(UA_Y);

    vX += dvX;
    vY += dvY;
    P += dP;

    double newerr = (GX*vX + GY*vY + dvX + dvY).squaredNorm()/(SIZE*SIZE);
    if(newerr > err) divergence = true;
  }

  if(divergence) std::cout<<"Instability encountered."<<std::endl;
}
