#include <noise/noise.h>

class Field{
public:
  //Grid-Size
  double dx = 1.0/(double)SIZE;
  double dy = 1.0/(double)SIZE;
  double dt = 0.005;

  //Parameters
  double viscosity = 0.001;     //[m^2/s]  //Future: Replace with Temperature Dependent?
  double density = 1.225;       //[kg/m^3]Future: Replace with Temperature Dependent?
  double P0 = 1E+5;             //Initial Pressure [Pa]

  //Volume Force
  glm::vec2 g = glm::vec2(0.0, 0.0);

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
  Eigen::SparseMatrix<double> B_MAT;
};

void Field::initialize(){
  //Boundary Vector
  /*
  boundary = Eigen::ArrayXd::Ones(SIZE*SIZE);
  boundary -= 0.9*shape::circle(glm::vec2(3.0*SIZE/4.0, SIZE/4.0), SIZE/10.0);
  boundary -= 0.9*shape::circle(glm::vec2(1.0*SIZE/4.0, SIZE/4.0), SIZE/10.0);
  boundary -= 0.9*shape::circle(glm::vec2(3.0*SIZE/4.0, 3.0*SIZE/4.0), SIZE/10.0);
  boundary -= 0.9*shape::circle(glm::vec2(1.0*SIZE/4.0, 3.0*SIZE/4.0), SIZE/10.0);

  boundary -= 0.9*shape::circle(glm::vec2(0.0*SIZE/4.0, 0.0*SIZE/4.0), SIZE/10.0);
  boundary -= 0.9*shape::circle(glm::vec2(2.0*SIZE/4.0, 0.0*SIZE/4.0), SIZE/10.0);
  boundary -= 0.9*shape::circle(glm::vec2(4.0*SIZE/4.0, 0.0*SIZE/4.0), SIZE/10.0);

  boundary -= 0.9*shape::circle(glm::vec2(0.0*SIZE/4.0, 2.0*SIZE/4.0), SIZE/10.0);
  boundary -= 0.9*shape::circle(glm::vec2(2.0*SIZE/4.0, 2.0*SIZE/4.0), SIZE/10.0);
  boundary -= 0.9*shape::circle(glm::vec2(4.0*SIZE/4.0, 2.0*SIZE/4.0), SIZE/10.0);

  boundary -= 0.9*shape::circle(glm::vec2(0.0*SIZE/4.0, 4.0*SIZE/4.0), SIZE/10.0);
  boundary -= 0.9*shape::circle(glm::vec2(2.0*SIZE/4.0, 4.0*SIZE/4.0), SIZE/10.0);
  boundary -= 0.9*shape::circle(glm::vec2(4.0*SIZE/4.0, 4.0*SIZE/4.0), SIZE/10.0);
  */

  boundary = Eigen::ArrayXd::Ones(SIZE*SIZE);
  boundary -= shape::circle(glm::vec2(SIZE/2.0, SIZE/2.0), SIZE/6.0);
  B_MAT = alg::sparseDiagonalize(boundary);

  //Initialize Pressure to some pressures pike in the center...
  P = P0*Eigen::ArrayXd::Ones(SIZE*SIZE);

  //Initialize Velocities to Zero
  vX = 0.0*Eigen::ArrayXd::Ones(SIZE*SIZE);
  vY = 1.0*Eigen::ArrayXd::Ones(SIZE*SIZE);
  vX = B_MAT*vX;
  vY = B_MAT*vY;

  PDE::initialize(dx, dy);
}

/*
  Solver Requirements:

    -> Solve the Pressure / Velocity Equations!!
    -> Solve the Transport Equation...

    First: Reduce the Algorithm Expressions in Complexity...
    i.e. Break the Algorithm down to its core elements...

    1. Define the Implicit Integration Matrix
    2.
*/

void Field::timestep(){

  Eigen::VectorXd E = Eigen::ArrayXd::Ones(SIZE*SIZE);

  //Construct Non-Linear Operator from Velocity Field and Other Contributions
  Eigen::SparseMatrix<double> vX_MAT = alg::sparseDiagonalize(vX);
  Eigen::SparseMatrix<double> vY_MAT = alg::sparseDiagonalize(vY);
  Eigen::SparseMatrix<double> E_MAT = -vX_MAT*PDE::XFLUX - vY_MAT*PDE::YFLUX + viscosity*(PDE::XDIFFUSION + PDE::YDIFFUSION);
  E_MAT = B_MAT*E_MAT;

  Eigen::VectorXd P_SOURCE_X = g.x*E - PDE::XFLUX*P;
  Eigen::VectorXd P_SOURCE_Y = g.y*E - PDE::YFLUX*P;
  P_SOURCE_X = B_MAT*P_SOURCE_X;
  P_SOURCE_Y = B_MAT*P_SOURCE_Y;

  //Implicit Solution
  Eigen::SparseMatrix<double> I_MAT_X = (alg::sparseIdentity() - dt*E_MAT);
  Eigen::SparseMatrix<double> I_MAT_Y = (alg::sparseIdentity() - dt*E_MAT);

  //All Elements!
  Eigen::VectorXd UA_X = E.cwiseQuotient(I_MAT_X*E);
  Eigen::VectorXd UA_Y = E.cwiseQuotient(I_MAT_Y*E);

  //Implicit Solution
  PDE::integrate<PDE::CN>(dt, vX, E_MAT, P_SOURCE_X);
  PDE::integrate<PDE::CN>(dt, vY, E_MAT, P_SOURCE_Y);
  PDE::integrate<PDE::EE>(dt, P, E_MAT);  //Guess the Pressure

  Eigen::VectorXd vXG;
  Eigen::VectorXd vYG;

  Eigen::VectorXd dP;
  Eigen::VectorXd dvX;
  Eigen::VectorXd dvY;

  double newerr = 1.0;
  double pCorr = 1.0;
  int maxiter = 100;
  int n = 0;

  while( (newerr > 1E-4 && pCorr > 1E-6) && !divergence && maxiter){
    n++;

    vX_MAT = alg::sparseDiagonalize(vX);
    vY_MAT = alg::sparseDiagonalize(vY);
    E_MAT = -vX_MAT*PDE::XFLUX - vY_MAT*PDE::YFLUX + viscosity*(PDE::XDIFFUSION + PDE::YDIFFUSION);
    E_MAT = B_MAT*E_MAT;

    //Compute the Intermediary Values
    P_SOURCE_X = g.x*E - PDE::XFLUX*P;
    P_SOURCE_Y = g.y*E - PDE::YFLUX*P;
    P_SOURCE_X = B_MAT*P_SOURCE_X;
    P_SOURCE_Y = B_MAT*P_SOURCE_Y;

    //E_MAT is updated already...
    I_MAT_X = (alg::sparseIdentity() - dt*E_MAT);
    I_MAT_Y = (alg::sparseIdentity() - dt*E_MAT);
    UA_X = E.cwiseQuotient(I_MAT_X*E);
    UA_Y = E.cwiseQuotient(I_MAT_Y*E);

    //Pressure Correction
    PDE::solver.compute(I_MAT_X);
    vXG.noalias() = PDE::solver.solve(vX+dt*P_SOURCE_X);
    PDE::solver.compute(I_MAT_Y);
    vYG.noalias() = PDE::solver.solve(vY+dt*P_SOURCE_Y);

    //
    PDE::solver.compute((PDE::GXF*alg::sparseDiagonalize(UA_X)*PDE::GXB + PDE::GYF*alg::sparseDiagonalize(UA_Y)*PDE::GYB)); //2D Laplace Operator
    dP = PDE::solver.solve((PDE::GXB*vXG + PDE::GYB*vYG));

    // We compute the velocity correction based on the pressure correction

    dvX = -(PDE::GXF*dP).cwiseProduct(UA_X);
    dvY = -(PDE::GYF*dP).cwiseProduct(UA_Y);

    // We correct our velocity guesses from the intermediary field

    vX = vXG + 0.9*B_MAT*dvX;
    vY = vYG + 0.9*B_MAT*dvY;

    // We correct our pressure guess from our previous guess
    P += B_MAT*dP;

    // Compute the Error (Divergence of Field) simply using our guess.
    newerr = (PDE::GXB*vXG + PDE::GYB*vYG + dvY + dvX).squaredNorm()/(SIZE*SIZE);
    pCorr = dP.squaredNorm()/(SIZE*SIZE);
    maxiter--;

    if(newerr > err) divergence = true;
  }
}
