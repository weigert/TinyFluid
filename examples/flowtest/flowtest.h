#include <noise/noise.h>

class Field{
public:
  //Topographical Information
  Eigen::VectorXd heightMap;
  double sealevel = 0.5;

  //Some Quantity
  Eigen::VectorXd conc;

  //Velocitites
  double u = 1.0;
  double v = 1.0;
  double D = 0.0001;

  double dx = 1.0/(double)SIZE;
  double dy = 1.0/(double)SIZE;

  //Parameters
  double viscosity = 0.0000148; //[m^2/s]  //Future: Replace with Temperature Dependent?
  double density = 1.225;       //[kg/m^3]Future: Replace with Temperature Dependent?
  double P0 = 1E+5; //Initial Pressure [Pa]
  double T0 = 300.0; //Initial Temperature

  //Cell Volume and Contact Area
  Eigen::VectorXd volume;

  //Volume Force

  //Integration Stuff
  void initialize();
  void timestep();
  int count = 0;


  Eigen::SparseMatrix<double> I;// = Eigen::MatrixXd::Identity(SIZE*SIZE,SIZE*SIZE);

  //Velocities (Actually)
  Eigen::VectorXd vX;
  Eigen::VectorXd vY;
  Eigen::VectorXd P;

  Eigen::SparseMatrix<double> vX_MAT;
  Eigen::SparseMatrix<double> vY_MAT;

  //Diffusivity Terms at Boundaries
  Eigen::SparseMatrix<double> dX_f;
  Eigen::SparseMatrix<double> dX_b;
  Eigen::SparseMatrix<double> dY_f;
  Eigen::SparseMatrix<double> dY_b;
  Eigen::SparseMatrix<double> J;
  Eigen::SparseMatrix<double> G;

  //Laplace Operators
  Eigen::SparseMatrix<double> LAPLACE;  //Full Laplace Operator in 2D
  Eigen::SparseMatrix<double> XFLUX;    //Integrated Surface Flux X, UNIFORM AREA!!
  Eigen::SparseMatrix<double> YFLUX;

  Eigen::SparseMatrix<double> GX;
  Eigen::SparseMatrix<double> GY;

  //Full Operator (Space-Discretization)!
  Eigen::SparseMatrix<double> MAT;
};



void Field::initialize(){

  //Initialize Heightmap
  noise::module::Perlin perlin;
  perlin.SetOctaveCount(20);
  perlin.SetFrequency(0.5);
  perlin.SetPersistence(0.8);

  heightMap = Eigen::ArrayXd::Zero(SIZE*SIZE);

  for(int i = 0; i < SIZE; i++){
    for(int j = 0; j < SIZE; j++){
      heightMap(am::index(i, j)) = (perlin.GetValue(i*(1.0/(double)SIZE), j*(1.0/(double)SIZE), 0.0) + 2.5) / 5.0;
    }
  }

  //Initialize Volume and Area Maps
  P = P0*Eigen::ArrayXd::Ones(SIZE*SIZE);

  //Initialize Velocities to Zero
  vX = u*Eigen::ArrayXd::Ones(SIZE*SIZE);
  vY = v*Eigen::ArrayXd::Ones(SIZE*SIZE);
  vX += u*am::flatGaussian();
  vY += v*am::flatGaussian();

  vX_MAT = am::sparseDiagonalize(vX);
  vY_MAT = am::sparseDiagonalize(vY);

  dX_f = am::getDiffX(1);
  dX_b = am::getDiffX(-1);
  dY_f = am::getDiffY(1);
  dY_b = am::getDiffY(-1);

  //Discrete Differential Operators (of some order)
  LAPLACE = am::laplaceX()/dx/dx + am::laplaceY()/dy/dy;
  GX = am::gradX()/dx;  //Gradient needs to be divided by the grid-spacing
  GY = am::gradY()/dy;  //-``-

  XFLUX = (am::getTrapezeX(1) - am::getTrapezeX(-1))/dx;
  YFLUX = (am::getTrapezeY(1) - am::getTrapezeY(-1))/dy;

  //Diffusion Operator! (Note: This is in FV Form and Includes the Area / Volume!!!)
  G = (dX_f+dX_b)/dx/dx + (dY_f+dY_b)/dy/dy;
}

void Field::timestep(){

  std::cout<<"Timestep: "<<count++<<std::endl;
  double dt = 0.1;

  /* This needs to be updated of course!! */
  vX_MAT = am::sparseDiagonalize(vX);
  vY_MAT = am::sparseDiagonalize(vY);

  //Include the Velocity Field here!
  Eigen::SparseMatrix<double> MAT = -vX_MAT*XFLUX -vY_MAT*YFLUX + viscosity*G;

  //Second Attempt: Choi et al. (1994) (NON-ITERATIVE)

  //Add General Source Term
  Eigen::VectorXd P_SOURCE_X = 0.0*Eigen::ArrayXd::Ones(SIZE*SIZE);
  Eigen::VectorXd P_SOURCE_Y = 0.0*Eigen::ArrayXd::Ones(SIZE*SIZE);

  //Add the Pressure Source
  P_SOURCE_X += -2.0*XFLUX*P;
  P_SOURCE_Y += -2.0*YFLUX*P;

  Eigen::VectorXd dP; //Pressure Correction!

  Eigen::VectorXd vX_Temp = vX;
  Eigen::VectorXd vY_Temp = vY;

  //Sparse Eigen-Solver
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(LAPLACE);                      //2D Laplace Operator

  //This section could be converged! (To get a better in-between estimate)
    //Crank-Nicholson Integration (With full Explicit Source)
    PDE::integrate<PDE::EE>(0.5*dt, vX_Temp, MAT, P_SOURCE_X);
    PDE::integrate<PDE::IE>(0.5*dt, vX_Temp, MAT);
    PDE::integrate<PDE::EE>(0.5*dt, vY_Temp, MAT, P_SOURCE_Y);
    PDE::integrate<PDE::IE>(0.5*dt, vY_Temp, MAT);

    //Pressure Correction
    dP = solver.solve((GX*vX_Temp + GY*vY_Temp)*2.0/dt); //Maybe this dt needs to go

  //Final Calculation (Updated Pressure Correction Source and Combine!)
  P_SOURCE_X = -XFLUX*dP;        //New Values from Pressure Correction
  P_SOURCE_Y = -YFLUX*dP;
  vX = vX_Temp + dt * 0.5 * P_SOURCE_X;     //Only Consider Pressure Correction!!
  vY = vY_Temp + dt * 0.5 * P_SOURCE_Y;

  //Fix the Pressure
  P += dP;

  //Finished!
}


  //Clear the Renderer
