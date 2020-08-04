#include <noise/noise.h>

class Field{
public:
  //Grid-Size
  double dx = 1.0/(double)SIZE;
  double dy = 1.0/(double)SIZE;
  double dt = 0.02;

  //Parameters
  double viscosity = 0.001;    //[m^2/s]  //Future: Replace with Temperature Dependent?
  double density = 1.225;       //[kg/m^3]Future: Replace with Temperature Dependent?
  double P0 = 1E+5;             //Initial Pressure [Pa]

  int SEED = 0;

  //Volume Force
  glm::vec2 g = glm::vec2(1.0, 1.0);
  glm::vec2 v0 = glm::vec2(1.0, 1.0); //Critical Velocity

  float sealevel;

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

  //Extra Fields
  Eigen::VectorXd humidity;
  Eigen::VectorXd temperature;
  double H0 = 0.01;  //Initial Humidity
  double T0 = 25.0; //Celsius

  //If these two guys are too large, the disturbances propagate like a wave
  //i.e. elliptically instead of characteristically
  //If their ratio is too different, then one propagates too fast or too slow
  //Meaning you get these huge eliminations or no spread of rain nucleation
  double diffusivity = 0.0004;
  double conductivity = 0.0008; //Important: Diffusivity Lower than Conductivity (Slightly)

  //Boundary Condition
  Eigen::SparseMatrix<double> B_MAT;
  Eigen::SparseMatrix<double> H_MAT;
  Eigen::SparseMatrix<double> V_MAT;
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
    height(i) = perlin.GetValue(pos.x/(double)SIZE, pos.y/(double)SIZE, SEED);
  }

  Eigen::VectorXd E = Eigen::ArrayXd::Ones(SIZE*SIZE);

  //Rescale the Height (First Multiply by Gaussian)
  height = height.cwiseProduct(shape::flatGaussian(glm::vec2(SIZE/2.0, SIZE/2.0), 4.0*SIZE));
  double hmin = height.minCoeff();
  double hmax = height.maxCoeff();
  height = (height - hmin*E)/(hmax - hmin);
  sealevel = 1.1*height.mean();

  //Compute the Boundary
  boundary = Eigen::ArrayXd::Zero(SIZE*SIZE);
  Eigen::VectorXd a = Eigen::ArrayXd::Zero(SIZE*SIZE);
  for(int i = 0; i < SIZE*SIZE; i++){
    //We can raise the boundary by altering this ratio...
    a(i) = (1.0-height(i))/(1.0-sealevel)*0.9;
    boundary(i) = (height(i) > sealevel)?a(i)*0.3:1.0;
  }

  //Square it for more emphasis on choking effect
  boundary = boundary.cwiseProduct(boundary);
  B_MAT = alg::sparseDiagonalize(boundary);
  H_MAT = alg::sparseDiagonalize(a);

  //Blower Intensity...
  Eigen::VectorXd b = Eigen::ArrayXd::Zero(SIZE*SIZE);
  b += shape::flatGaussian(glm::vec2(0.0, SIZE/2.0), SIZE/2.0);
  b += shape::flatGaussian(glm::vec2(SIZE, SIZE/2.0), SIZE/2.0);
  V_MAT = alg::sparseDiagonalize(b); //VMAT Specifies exactly where volume forces apply.

  //Set humidity
  humidity = 0.0*E;
  temperature = T0*E;

  //Initialize Pressure to some pressures pike in the center...
  P = P0*E;
  //P += 0.0*(E-B_MAT*E); //Pressure Proportional to distance from sealevel?

  //Initialize Velocities to Zero
  vX = Eigen::ArrayXd::Zero(SIZE*SIZE);
  vY = Eigen::ArrayXd::Zero(SIZE*SIZE);
  vX = B_MAT*vX;
  vY = B_MAT*vY;

  PDE::initialize(dx, dy);
}

bool turned = false;

Eigen::VectorXd E = Eigen::ArrayXd::Ones(SIZE*SIZE);
Eigen::SparseMatrix<double> E_MAT;
Eigen::SparseMatrix<double> P_MAT;

Eigen::VectorXd P_SOURCE_X;
Eigen::VectorXd P_SOURCE_Y;

Eigen::VectorXd vXG;
Eigen::VectorXd vYG;

Eigen::VectorXd dP;
Eigen::VectorXd dvX;
Eigen::VectorXd dvY;

Eigen::VectorXd A_X;
Eigen::VectorXd A_Y;

Eigen::VectorXd UA_X;
Eigen::VectorXd UA_Y;

Eigen::SparseMatrix<double> TRANSPORTMAT_H;
Eigen::SparseMatrix<double> TRANSPORTMAT_T;

Eigen::SparseMatrix<double> I_MAT_X;
Eigen::SparseMatrix<double> I_MAT_Y;

Eigen::VectorXd HSOURCE;
Eigen::VectorXd TSOURCE;

void Field::timestep(){
  count++;
  int time = 1000;
  if(count > time/2 && !turned) { //Everytime we are half way through a cycle...
    v0 = glm::vec2(-v0.y, v0.x); //Wind Shifts by 90 Degrees!
    turned = true;
  }
  if(count > time){
    turned = false;
    count = 0;
  }

  //Construct Non-Linear Operator from Velocity Field and Other Contributions
  E_MAT = -1.0*(vX.asDiagonal()*PDE::XFLUX + vY.asDiagonal()*PDE::YFLUX) + viscosity*(PDE::XDIFFUSION + PDE::YDIFFUSION);
  E_MAT = B_MAT*E_MAT;
  P_MAT = -PDE::XFLUX -PDE::YFLUX;
  P_MAT = B_MAT*P_MAT;

  P_SOURCE_X = g.x*V_MAT*(E*v0.x-vX) - PDE::XFLUX*P;
  P_SOURCE_Y = g.y*V_MAT*(E*v0.y-vY) - PDE::YFLUX*P;
  P_SOURCE_X = B_MAT*P_SOURCE_X;
  P_SOURCE_Y = B_MAT*P_SOURCE_Y;

  //Implicit Solution
  I_MAT_X = (alg::sparseIdentity() - dt*E_MAT);
  I_MAT_Y = (alg::sparseIdentity() - dt*E_MAT);

  //Helper Values (Sum over Rows)
  A_X = I_MAT_X*E;
  A_Y = I_MAT_Y*E;

  //All Elements!
  UA_X = E.cwiseQuotient(A_X);
  UA_Y = E.cwiseQuotient(A_Y);

  //Implicit Solution
  PDE::integrate<PDE::CN>(dt, vX, E_MAT, P_SOURCE_X);
  PDE::integrate<PDE::CN>(dt, vY, E_MAT, P_SOURCE_Y);

//  TRANSPORTMAT_H = H_MAT*(-1.0*(vX.asDiagonal()*PDE::XFLUX + vY.asDiagonal()*PDE::YFLUX) + diffusivity*(PDE::XDIFFUSION + PDE::YDIFFUSION));
//  TRANSPORTMAT_T = H_MAT*(-1.0*(vX.asDiagonal()*PDE::XFLUX + vY.asDiagonal()*PDE::YFLUX) + conductivity*(PDE::XDIFFUSION + PDE::YDIFFUSION));
//  HSOURCE = source::HSOURCE(boundary, height, sealevel, humidity, P, temperature);
//  TSOURCE = source::TSOURCE(boundary, height, sealevel, humidity, P, temperature, vX, vY);
//  PDE::integrate<PDE::CN>(dt, humidity, TRANSPORTMAT_H, HSOURCE);
//  PDE::integrate<PDE::CN>(dt, temperature, TRANSPORTMAT_T, TSOURCE);

  double newerr = 1.0;
  double pCorr = 1.0;
  int maxiter = 100;
  int n = 0;

  while(newerr > 1E-4 && pCorr > 1E-6 && !divergence && maxiter){
    n++;

    E_MAT = -1.0*(vX.asDiagonal()*PDE::XFLUX + vY.asDiagonal()*PDE::YFLUX) + viscosity*(PDE::XDIFFUSION + PDE::YDIFFUSION);
    E_MAT = B_MAT*E_MAT;

    //Compute the Intermediary Values
    P_SOURCE_X = g.x*V_MAT*(E*v0.x-vX) - PDE::XFLUX*P;
    P_SOURCE_Y = g.y*V_MAT*(E*v0.y-vY) - PDE::YFLUX*P;
    P_SOURCE_X = B_MAT*P_SOURCE_X;
    P_SOURCE_Y = B_MAT*P_SOURCE_Y;

    //E_MAT is updated already...
    I_MAT_X = (alg::sparseIdentity() - dt*E_MAT);
    I_MAT_Y = (alg::sparseIdentity() - dt*E_MAT);
    A_X = I_MAT_X*E;
    A_Y = I_MAT_Y*E;
    UA_X = E.cwiseQuotient(A_X);
    UA_Y = E.cwiseQuotient(A_Y);

    //Pressure Correction

    PDE::solver.compute(I_MAT_X);
    vXG.noalias() = PDE::solver.solve(vX+dt*P_SOURCE_X);
    PDE::solver.compute(I_MAT_Y);
    vYG.noalias() = PDE::solver.solve(vY+dt*P_SOURCE_Y);

    //Pressure Correction
    PDE::solver.factorize((PDE::GXF*alg::sparseDiagonalize(UA_X)*PDE::GXB+ PDE::GYF*alg::sparseDiagonalize(UA_Y)*PDE::GYB)); //2D Laplace Operator
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

    if(newerr > err) divergence = true;
    maxiter--;

  //  TRANSPORTMAT_H = H_MAT*(-1.0*(vX.asDiagonal()*PDE::XFLUX + vY.asDiagonal()*PDE::YFLUX) + diffusivity*(PDE::XDIFFUSION + PDE::YDIFFUSION));
  //  TRANSPORTMAT_T = H_MAT*(-1.0*(vX.asDiagonal()*PDE::XFLUX + vY.asDiagonal()*PDE::YFLUX) + conductivity*(PDE::XDIFFUSION + PDE::YDIFFUSION));
  //  HSOURCE = source::HSOURCE(boundary, height, sealevel, humidity, P, temperature);
  //  TSOURCE = source::TSOURCE(boundary, height, sealevel, humidity, P, temperature, vX, vY);
  //  PDE::integrate<PDE::CN>(dt, humidity, TRANSPORTMAT_H, HSOURCE);
  //  PDE::integrate<PDE::CN>(dt, temperature, TRANSPORTMAT_T, TSOURCE);
  }

  //Transport...
}
