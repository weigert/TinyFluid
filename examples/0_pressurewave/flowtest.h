#include <noise/noise.h>

class Field{
public:
  //Grid-Size
  double dx = 1.0/(double)SIZE;
  double dy = 1.0/(double)SIZE;
  double dt = 0.0001;

  //Parameters
  double viscosity = 0.005;    //[m^2/s]  //Future: Replace with Temperature Dependent?
  double density = 1.225;       //[kg/m^3]Future: Replace with Temperature Dependent?
  double P0 = 1E+5;             //Initial Pressure [Pa]

  void initialize();
  void timestep();

  //Velocities (Actually)
  Eigen::VectorXd vX;
  Eigen::VectorXd vY;
  Eigen::VectorXd P;

  Eigen::VectorXd boundary;
  Eigen::VectorXd B;

  glm::vec2 g;
  glm::vec2 v0;
};

void Field::initialize(){

  //Initialize Fields
  vX = Eigen::ArrayXd::Zero(SIZE*SIZE);
  vY = Eigen::ArrayXd::Zero(SIZE*SIZE);
  P = P0*Eigen::ArrayXd::Ones(SIZE*SIZE);

  //Pressure Wave
  P += (50.0*shape::flatGaussian(glm::vec2(SIZE/2.0, SIZE/2.0), SIZE/3.0)).matrix();

  //Boundary Conditions
  boundary = Eigen::ArrayXd::Ones(SIZE*SIZE);

  //Blower
  B = Eigen::ArrayXd::Zero(SIZE*SIZE);
  g = glm::vec2(0);
  v0 = glm::vec2(0);

  PDE::initialize(dx, dy);

}

void Field::timestep(){

  //Solve the Navier Stokes Equations
  PDE::navierstokes(dt, viscosity, vX, vY, P, boundary, g, v0, B);

}
