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
  double viscosity = 0.000148; //[m^2/s]  //Future: Replace with Temperature Dependent?
  double density = 1.225;       //[kg/m^3]Future: Replace with Temperature Dependent?
  double P0 = 1E+3; //Initial Pressure [Pa]
  double T0 = 300.0; //Initial Temperature

  //Cell Volume and Contact Area
  Eigen::VectorXd volume;

  //Integration Stuff
  void initialize();
  void timestep();
  int count = 0;

  //Velocities (Actually)
  Eigen::VectorXd vX;
  Eigen::VectorXd vY;
  Eigen::VectorXd P;

  //Laplace Operator
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  Eigen::SparseMatrix<double> LAPLACE;  //Full Laplace Operator in 2D
  Eigen::MatrixXd LAPLACE_DENSE;

  //Gradient
  Eigen::SparseMatrix<double> GX;       //Gradient Operator X-Element
  Eigen::SparseMatrix<double> GY;       //Gradient Operator Y-Element

  Eigen::SparseMatrix<double> XFLUX;         //Integrated Surface Flux X, UNIFORM AREA!!
  Eigen::SparseMatrix<double> YFLUX;
  Eigen::SparseMatrix<double> XDIFFUSION;    //Integrated Surface Flux X, UNIFORM AREA!!
  Eigen::SparseMatrix<double> YDIFFUSION;
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

  //Initialize Pressure to Something Non-Obvious
  P = P0*Eigen::ArrayXd::Ones(SIZE*SIZE);
  P += 100*am::flatGaussian();

  //Initialize Velocities to Zero
  vX = Eigen::ArrayXd::Zero(SIZE*SIZE);
  vY = Eigen::ArrayXd::Zero(SIZE*SIZE);

  //Discrete Differential Operators (of some order)
  GX = am::gradX()/dx;  //Gradient needs to be divided by the grid-spacing
  GY = am::gradY()/dy;  //-``-

  LAPLACE = am::laplaceX()/dx/dx + am::laplaceY()/dy/dy;
  LAPLACE_DENSE = LAPLACE;

  XFLUX = (am::getTrapezeX(1) - am::getTrapezeX(-1))/dx;
  YFLUX = (am::getTrapezeY(1) - am::getTrapezeY(-1))/dy;

  //Diffusion Operator! (Note: This is in FV Form and Includes the Area / Volume!!!)
  XDIFFUSION = (am::getDiffX(1)+am::getDiffX(-1))/dx/dx;
  YDIFFUSION = (am::getDiffY(1)+am::getDiffY(-1))/dy/dy;
}

void Field::timestep(){
  std::cout<<"Timestep: "<<count++<<std::endl;
  double dt = 0.005;

  /* SIMPLEC ALGORITHM */

  //Get the Original Values for
  Eigen::VectorXd vX_Temp = vX;
  Eigen::VectorXd vY_Temp = vY;
  Eigen::VectorXd P_Temp = P;

  Eigen::VectorXd zero = Eigen::ArrayXd::Zero(SIZE*SIZE);

  // I am not sure if these operators need to be redfined for the new vX_Temp or not!!!

  //Construct Non-Linear Operator from Velocity Field and Other Contributions
  Eigen::SparseMatrix<double> vX_MAT = am::sparseDiagonalize(vX);
  Eigen::SparseMatrix<double> vY_MAT = am::sparseDiagonalize(vY);

  //Implicit and Explicitly Evaluated Matrices
  Eigen::SparseMatrix<double> I_MAT_X = -vX_MAT*XFLUX + viscosity*XDIFFUSION;
  Eigen::SparseMatrix<double> I_MAT_Y = -vY_MAT*YFLUX + viscosity*YDIFFUSION;
  Eigen::SparseMatrix<double> I_MAT = I_MAT_X+I_MAT_Y;
  Eigen::SparseMatrix<double> E_MAT = am::sparseDiagonalize(zero);

  // MAT could theoretically be split into an implicit and explicit portion

  //Construct the weight vector??? ONLY FROM THE IMPLICIT MATRIX PORTION!!
  Eigen::VectorXd E = Eigen::ArrayXd::Ones(SIZE*SIZE);

  //All Elements! (Summed for Each Individual Guy)
  Eigen::VectorXd A_X = ((am::sparseIdentity()-dt*I_MAT_X)*E);
  Eigen::VectorXd A_Y = ((am::sparseIdentity()-dt*I_MAT_Y)*E);

  //All Elements!
  Eigen::VectorXd UA_X = E.cwiseQuotient(A_X);
  Eigen::VectorXd UA_Y = E.cwiseQuotient(A_Y);

  //The solver needs to be set up with a specific matrix... this one!
  //solver.compute(GX*am::sparseDiagonalize(UA_X)*GX + GY*am::sparseDiagonalize(UA_Y)*GY);                      //2D Laplace Operator
  solver.compute(am::sparseDiagonalize(UA_X)*am::laplaceX()/dx/dx + am::sparseDiagonalize(UA_Y)*am::laplaceY()/dy/dy);                      //2D Laplace Operator

  Eigen::VectorXd dP;
  Eigen::VectorXd dvX;
  Eigen::VectorXd dvY;

  Eigen::VectorXd P_SOURCE_X;
  Eigen::VectorXd P_SOURCE_Y;

  //Iterate!!
  for(int i = 0; i < 10; i++){

    //Pressure Source-Term from Previous Iteration
     P_SOURCE_X = -XFLUX*P;
     P_SOURCE_Y = -YFLUX*P;

    //Solution of Momentum Equations with Underrelaxation! Including the Pressure Source.

    //Notes: All matrices (Implicit + Explicit) should be constructed from the original vX of the previous "inner iteration".
    //Only the vector they then act upon is the outer iteration portion. Also the pressure source is from the outer iteration.
    //Outer Iteration = Convergence Step, Inner Iteration = Timestep

    //To use a point implicit method you need to define the explicit matrix also with the original vX
    //But then let it act on the vX of the previous iteration!!!! So define it outside this loop.

    //That is why MAT is not changed inside the inner iteration, but the source term here is!!
    //PDE::integrate<PDE::CN>(dt, vX, I_MAT, P_SOURCE_X);
    //PDE::integrate<PDE::CN>(dt, vY, I_MAT, P_SOURCE_Y);

    PDE::integrate<PDE::EE>(dt, vX, E_MAT, P_SOURCE_X);
    PDE::integrate<PDE::IE>(dt, vX, I_MAT);
    PDE::integrate<PDE::EE>(dt, vY, E_MAT, P_SOURCE_Y);
    PDE::integrate<PDE::IE>(dt, vY, I_MAT);

    //Pressure and Velocity Correction (Intermediate)
    dP = solver.solve((GX*vX + GY*vY));
    dvX = -(GX*dP).cwiseProduct(UA_X);
    dvY = -(GY*dP).cwiseProduct(UA_Y);

    //Correct the Temporary Values!
    //Underrelaxation is necessary for this thing to work properly.
    //This raises convergence speed. (Because of numerical overshooting I guess)
    //If you don't underrelax enough (i.e. alpha = 0.2) then it will give very steep increases
    //0.5 is much better... we can reach a divergence free field after only a few sets of 20 iterations...
    //If our field is not divergence free then we are generally unstable...
    //0.7 is even better.
    //
    //1.0 doesn't work at all!
    //A good underrelaxation value for the speed is ~0.8

    vX += 0.7*dvX;
    vY += 0.7*dvY;
    P += dP;

    //Convergence Criterion!
    //std::cout << (GX*vX + GY*vY + dvX + dvY).squaredNorm() << std::endl;
    //std::cout<<dP<<std::endl;
  }
}
