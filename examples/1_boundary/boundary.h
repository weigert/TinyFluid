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
  glm::vec2 g = glm::vec2(1.0, 0.0);  //Volume Force
  glm::vec2 v0 = glm::vec2(1.0, 0.0); //Terminal Velocity

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

  Eigen::VectorXd B; //Boundary
  Eigen::VectorXd F; //Volume Force
};

void Field::initialize(){

  //Boundary Condition
  B = Eigen::ArrayXd::Ones(SIZE*SIZE);
  B -= shape::circle(glm::vec2(SIZE/2.0, SIZE/2.0), SIZE/6.0);

  //Pressure and Velocity Fields
  P = P0*Eigen::ArrayXd::Ones(SIZE*SIZE);
  vX = 1.0*Eigen::ArrayXd::Ones(SIZE*SIZE);
  vY = 0.0*Eigen::ArrayXd::Ones(SIZE*SIZE);
  vX = B.cwiseProduct(vX);
  vY = B.cwiseProduct(vY);

  //Volume Forces
  F = Eigen::ArrayXd::Ones(SIZE*SIZE);

  //Setup Matrices
  PDE::initialize(dx, dy);
}

void Field::timestep(){

  PDE::navierstokes(dt, viscosity, vX, vY, P, B, g, v0, F);

}
