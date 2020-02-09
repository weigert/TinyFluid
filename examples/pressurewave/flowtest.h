#include <noise/noise.h>

class Field{
public:
  //Grid-Size
  double dx = 1.0/(double)SIZE;
  double dy = 1.0/(double)SIZE;
  double dt = 0.0005;

  //Parameters
  double viscosity = 0.000148;  //[m^2/s]  //Future: Replace with Temperature Dependent?
  double density = 1.225;       //[kg/m^3]Future: Replace with Temperature Dependent?
  double P0 = 1E+5;             //Initial Pressure [Pa]

  //Volume Force
  glm::vec2 g = glm::vec2(0.0);

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

  Eigen::VectorXd boundary;

  //Laplace Operator
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;

  //Gradient
  Eigen::SparseMatrix<double> GXF;       //Gradient Operator X-Element
  Eigen::SparseMatrix<double> GYF;       //Gradient Operator Y-Element
  Eigen::SparseMatrix<double> GXB;       //Gradient Operator X-Element
  Eigen::SparseMatrix<double> GYB;       //Gradient Operator Y-Element

  //Spatial Discretization Operator
  Eigen::SparseMatrix<double> XFLUX;         //Integrated Surface Flux X, UNIFORM AREA!!
  Eigen::SparseMatrix<double> YFLUX;
  Eigen::SparseMatrix<double> XDIFFUSION;    //Integrated Surface Flux X, UNIFORM AREA!!
  Eigen::SparseMatrix<double> YDIFFUSION;

  Eigen::SparseMatrix<double> B_MAT;

};

void Field::initialize(){
  //Initialize Pressure to some pressures pike in the center...
  P = P0*Eigen::ArrayXd::Ones(SIZE*SIZE);
  //P += 100*shape::flatGaussian(glm::vec2(SIZE/2.0, SIZE/2.0), SIZE/2.0);

  //Double Pressure Wave
  P += 100*shape::flatGaussian(glm::vec2(SIZE/3.0, SIZE/2.0), SIZE/3.0);
  P += 100*shape::flatGaussian(glm::vec2(2.0*SIZE/3.0, SIZE/2.0), SIZE/3.0);

  //Initialize Velocities to Zero
  vX = Eigen::ArrayXd::Zero(SIZE*SIZE);
  vY = Eigen::ArrayXd::Zero(SIZE*SIZE);

  boundary = Eigen::ArrayXd::Ones(SIZE*SIZE);
  B_MAT = alg::sparseDiagonalize(boundary);

  //Discrete Differential Operators (of some order)
  GXF = space::CFD(glm::vec2(1, 0), 1)/dx;  //Gradient needs to be divided by the grid-spacing
  GYF = space::CFD(glm::vec2(0, 1), 1)/dy;  //-``-
  GXB = space::SFD(glm::vec2(1, 0), 1)/dx;  //Gradient needs to be divided by the grid-spacing
  GYB = space::SFD(glm::vec2(0, 1), 1)/dy;  //-``-

  solver.compute(GXF*GXB + GYF*GYB);
  solver.setTolerance(0.000001);

  //Surface Integrators (i.e. convective fluxes over surface)
  XFLUX = space::FV_FLUX(glm::vec2(1, 0))/dx;
  YFLUX = space::FV_FLUX(glm::vec2(0, 1))/dy;

  //Diffusion Operator! (Note: This is in FV Form and Includes the Area / Volume!!!)
  XDIFFUSION = space::FV_DIFFUSION(glm::vec2(1, 0))/dx/dx;
  YDIFFUSION = space::FV_DIFFUSION(glm::vec2(0, 1))/dy/dy;
}

void Field::timestep(){

  /* SIMPLEC ALGORITHM
    -> Purely Explicit Here...
  */

  /*
    Approximate the Solution at the Next Time Point:
      -> Without the corrected pressure field!
      -> We guess the pressure field remains unchanged for now.
  */

  //Construct Non-Linear Operator from Velocity Field and Other Contributions
  Eigen::SparseMatrix<double> vX_MAT = alg::sparseDiagonalize(vX);
  Eigen::SparseMatrix<double> vY_MAT = alg::sparseDiagonalize(vY);
  Eigen::SparseMatrix<double> E_MAT = -vX_MAT*XFLUX - vY_MAT*YFLUX + viscosity*(XDIFFUSION + YDIFFUSION);

  Eigen::VectorXd P_SOURCE_X;
  Eigen::VectorXd P_SOURCE_Y;
  P_SOURCE_X = g.x*Eigen::ArrayXd::Ones(SIZE*SIZE);
  P_SOURCE_Y = g.y*Eigen::ArrayXd::Ones(SIZE*SIZE);
  P_SOURCE_X += -XFLUX*P;
  P_SOURCE_Y += -YFLUX*P;
  vX = vX + dt*B_MAT*(E_MAT*vX + P_SOURCE_X);
  vY = vY + dt*B_MAT*(E_MAT*vY + P_SOURCE_Y);

  Eigen::VectorXd vXG;
  Eigen::VectorXd vYG;

  /*
    Fix the pressure and velocity guesses.
  */

  //PDE::solveNavierStokes<PDE::SIMPLEC>(vX, vY, P, B_MAT, GXF, GXB, GYF, GYB, 1E-4, 100, 0.3, 0.9);

  Eigen::VectorXd dP;
  Eigen::VectorXd dvX;
  Eigen::VectorXd dvY;

  double newerr = 1.0;
  double pCorr = 1.0;
  int maxiter = 10000;
  int n = 0;

  while(newerr > 1E-4 && pCorr > 1E-4 && !divergence && maxiter){
    n++;

    vX_MAT = alg::sparseDiagonalize(vX);
    vY_MAT = alg::sparseDiagonalize(vY);
    E_MAT = -vX_MAT*XFLUX - vY_MAT*YFLUX + viscosity*(XDIFFUSION + YDIFFUSION);

    //Compute the Intermediary Values
    P_SOURCE_X = g.x*Eigen::ArrayXd::Ones(SIZE*SIZE);
    P_SOURCE_Y = g.y*Eigen::ArrayXd::Ones(SIZE*SIZE);
    P_SOURCE_X += -XFLUX*P;
    P_SOURCE_Y += -YFLUX*P;
    vXG = vX + dt*B_MAT*(E_MAT*vX + P_SOURCE_X);
    vYG = vY + dt*B_MAT*(E_MAT*vY + P_SOURCE_Y);

    //Pressure Correction

    dP = solver.solve((GXB*vXG + GYB*vYG));

    // We compute the velocity correction based on the pressure correction

    dvX = -(GXF*dP);
    dvY = -(GYF*dP);

    // We correct our velocity guesses from the intermediary field

    vX = vXG + 0.9*B_MAT*dvX;
    vY = vYG + 0.9*B_MAT*dvY;

    // We correct our pressure guess from our previous guess
    P += B_MAT*dP;

    // Compute the Error (Divergence of Field) simply using our guess.
    newerr = (GXB*vXG + GYB*vYG + dvY + dvX).squaredNorm()/(SIZE*SIZE);
    pCorr = dP.squaredNorm()/(SIZE*SIZE);
    maxiter--;

    if(newerr > err) divergence = true;
  }
}
