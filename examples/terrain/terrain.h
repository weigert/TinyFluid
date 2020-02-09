#include <noise/noise.h>

class Field{
public:
  //Grid-Size
  double dx = 1.0/(double)SIZE;
  double dy = 1.0/(double)SIZE;
  double dt = 0.01;

  //Parameters
  double viscosity = 0.000148;  //[m^2/s]  //Future: Replace with Temperature Dependent?
  double density = 1.225;       //[kg/m^3]Future: Replace with Temperature Dependent?
  double P0 = 1E+5;             //Initial Pressure [Pa]

  //Volume Force
  glm::vec2 g = glm::vec2(0.5, 0.5);

  float sealevel = 0.3;

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

  Eigen::VectorXd height;
  Eigen::VectorXd boundary;
  Eigen::VectorXd concentration;
  double c0 = 100;
  double diffusivity = 0.001;

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
  //Height Map Generation
  noise::module::Perlin perlin;
  perlin.SetOctaveCount(10);
  perlin.SetFrequency(6);
  perlin.SetPersistence(0.5);
  height = Eigen::ArrayXd::Zero(SIZE*SIZE);
  for(int i = 0; i < SIZE*SIZE; i++){
    glm::vec2 pos = alg::pos(i);
    height(i) = perlin.GetValue(pos.x/(double)SIZE, pos.y/(double)SIZE, 1.0);
  }

  Eigen::VectorXd E = Eigen::ArrayXd::Ones(SIZE*SIZE);

  //Rescale the Height
  double hmin = height.minCoeff();
  double hmax = height.maxCoeff();
  height = (height - hmin*E)/(hmax - hmin);

  //Compute the Boundary
  boundary = Eigen::ArrayXd::Zero(SIZE*SIZE);
  for(int i = 0; i < SIZE*SIZE; i++){
    boundary(i) = (height(i) > sealevel)?(1.0-height(i))/(1.0-sealevel):1.0;
  }

  B_MAT = alg::sparseDiagonalize(boundary);

  //Set Concentration
  concentration = c0*B_MAT*E;

  //Initialize Pressure to some pressures pike in the center...
  P = P0*Eigen::ArrayXd::Ones(SIZE*SIZE);

  //Initialize Velocities to Zero
  vX = 0.5*Eigen::ArrayXd::Ones(SIZE*SIZE);
  vY = 0.5*Eigen::ArrayXd::Ones(SIZE*SIZE);
  vX = B_MAT*vX;
  vY = B_MAT*vY;

  //Discrete Differential Operators (of some order)
  GXF = space::CFD(glm::vec2(1, 0), 1)/dx;  //Gradient needs to be divided by the grid-spacing
  GYF = space::CFD(glm::vec2(0, 1), 1)/dy;  //-``-
  GXB = space::SFD(glm::vec2(1, 0), 1)/dx;  //Gradient needs to be divided by the grid-spacing
  GYB = space::SFD(glm::vec2(0, 1), 1)/dy;  //-``-

  //Discrete Differential Operator (Laplace)
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

  Eigen::VectorXd E = Eigen::ArrayXd::Ones(SIZE*SIZE);


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

  //Implicit Solution
  Eigen::SparseMatrix<double> I_MAT_X = (alg::sparseIdentity() - dt*B_MAT*E_MAT);
  Eigen::SparseMatrix<double> I_MAT_Y = (alg::sparseIdentity() - dt*B_MAT*E_MAT);

  //Helper Values (Sum over Rows)
  Eigen::VectorXd A_X = I_MAT_X*E;
  Eigen::VectorXd A_Y = I_MAT_Y*E;

  //All Elements!
  Eigen::VectorXd UA_X = E.cwiseQuotient(A_X);
  Eigen::VectorXd UA_Y = E.cwiseQuotient(A_Y);

  //Implicit Solution
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver2;
  solver2.compute(I_MAT_X);
  vX = solver2.solve(vX+dt*B_MAT*P_SOURCE_X);
  solver2.compute(I_MAT_Y);
  vY = solver2.solve(vY+dt*B_MAT*P_SOURCE_Y);

/*
  vX = vX + dt*B_MAT*(E_MAT*vX + P_SOURCE_X);
  vY = vY + dt*B_MAT*(E_MAT*vY + P_SOURCE_Y);
*/


  Eigen::SparseMatrix<double> TRANSPORTMAT = B_MAT*(-vX_MAT*XFLUX - vY_MAT*YFLUX + diffusivity*(XDIFFUSION + YDIFFUSION));
  PDE::integrate<PDE::CN>(dt, concentration, TRANSPORTMAT);

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
  int maxiter = 100;
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

    //E_MAT is updated already...
    I_MAT_X = (alg::sparseIdentity() - dt*B_MAT*E_MAT);
    I_MAT_Y = (alg::sparseIdentity() - dt*B_MAT*E_MAT);
    A_X = I_MAT_X*E;
    A_Y = I_MAT_Y*E;
    UA_X = E.cwiseQuotient(A_X);
    UA_Y = E.cwiseQuotient(A_Y);

    //Pressure Correction
    solver2.compute(I_MAT_X);
    vXG = solver2.solve(vX+dt*B_MAT*P_SOURCE_X);
    solver2.compute(I_MAT_Y);
    vYG = solver2.solve(vY+dt*B_MAT*P_SOURCE_Y);

    /*
    vXG = vX;// + dt*B_MAT*(E_MAT*vX + P_SOURCE_X);
    vYG = vY;// + dt*B_MAT*(E_MAT*vY + P_SOURCE_Y);
*/
    //Pressure Correction
    solver.factorize((GXF*alg::sparseDiagonalize(UA_X)*GXB+ GYF*alg::sparseDiagonalize(UA_Y)*GYB)); //2D Laplace Operator
    dP = solver.solve((GXB*vXG + GYB*vYG));

    // We compute the velocity correction based on the pressure correction

    dvX = -(GXF*dP).cwiseProduct(UA_X);
    dvY = -(GYF*dP).cwiseProduct(UA_Y);

    // We correct our velocity guesses from the intermediary field

    vX = vXG + 0.9*B_MAT*dvX;
    vY = vYG + 0.9*B_MAT*dvY;

    // We correct our pressure guess from our previous guess
    P += B_MAT*dP;

    // Compute the Error (Divergence of Field) simply using our guess.
    newerr = (GXB*vXG + GYB*vYG + dvY + dvX).squaredNorm()/(SIZE*SIZE);
    pCorr = dP.squaredNorm()/(SIZE*SIZE);

    if(newerr > err) divergence = true;
    maxiter--;

    TRANSPORTMAT = B_MAT*(-vX_MAT*XFLUX - vY_MAT*YFLUX + diffusivity*(XDIFFUSION + YDIFFUSION));
    PDE::integrate<PDE::CN>(dt, concentration, TRANSPORTMAT);
  }

  //Transport...

}

/*
void Field::timestep(){

  Eigen::VectorXd E = Eigen::ArrayXd::Ones(SIZE*SIZE);

  //Construct Non-Linear Operator from Velocity Field and Other Contributions
  Eigen::SparseMatrix<double> vX_MAT = alg::sparseDiagonalize(vX);
  Eigen::SparseMatrix<double> vY_MAT = alg::sparseDiagonalize(vY);

  //Implicit and Explicitly Evaluated Operators...
  Eigen::SparseMatrix<double> E_MAT_X = -vX_MAT*XFLUX + viscosity*XDIFFUSION;
  Eigen::SparseMatrix<double> E_MAT_Y = -vY_MAT*YFLUX + viscosity*YDIFFUSION;

  Eigen::SparseMatrix<double> E_MAT = (E_MAT_X+E_MAT_Y);

  Eigen::SparseMatrix<double> I_MAT_X = 0.0*E_MAT_X;
  Eigen::SparseMatrix<double> I_MAT_Y = 0.0*E_MAT_Y;
  Eigen::SparseMatrix<double> I_MAT = 0.0*E_MAT;

  //All Elements! (Summed for Each Individual Guy)
  Eigen::VectorXd A_X = ((alg::sparseIdentity()-dt*I_MAT_X)*E);
  Eigen::VectorXd A_Y = ((alg::sparseIdentity()-dt*I_MAT_Y)*E);

  //All Elements!
  Eigen::VectorXd UA_X = E.cwiseQuotient(A_X);
  Eigen::VectorXd UA_Y = E.cwiseQuotient(A_Y);

  //Refactorize the Solver Matrix
  solver.factorize((GXF*alg::sparseDiagonalize(UA_X)*GXB+ GYF*alg::sparseDiagonalize(UA_Y)*GYB)); //2D Laplace Operator


  Eigen::VectorXd P_SOURCE_X;
  Eigen::VectorXd P_SOURCE_Y;
  P_SOURCE_X = g.x*Eigen::ArrayXd::Ones(SIZE*SIZE);
  P_SOURCE_Y = g.y*Eigen::ArrayXd::Ones(SIZE*SIZE);
  P_SOURCE_X += -XFLUX*P;
  P_SOURCE_Y += -YFLUX*P;
  vX = vX + dt*B_MAT*(E_MAT*vX + P_SOURCE_X);
  vY = vY + dt*B_MAT*(E_MAT*vY + P_SOURCE_Y);

  Eigen::VectorXd dP;
  Eigen::VectorXd dvX;
  Eigen::VectorXd dvY;

  double newerr = 1.0;
  int maxiter = 100;

  while(newerr > 1E-4 && !divergence && maxiter){

    //Pressure Correction

    dP = solver.solve((GXB*vX + GYB*vY));

    // We compute the velocity correction based on the pressure correction

    dvX = -(GXF*dP).cwiseProduct(UA_X);
    dvY = -(GYF*dP).cwiseProduct(UA_Y);

    // We correct our velocity guesses from the intermediary field

    vX += 0.9*B_MAT*dvX;
    vY += 0.9*B_MAT*dvY;

    // We correct our pressure guess from our previous guess

    P += 0.3*B_MAT*dP;

    // Compute the Error (Divergence of Field) simply using our guess.
    newerr = (GXB*vX + GYB*vY).squaredNorm()/(SIZE*SIZE);
    maxiter--;

    if(newerr > err) divergence = true;
  }
  if(!maxiter) std::cout<<"Max iterations surpassed."<<std::endl;
  if(divergence) std::cout<<"Instability encountered."<<std::endl;
}
*/
