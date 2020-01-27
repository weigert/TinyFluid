/*

  Time Integrators

*/

namespace PDE{

  enum Integrator{
    EE, //Explicit Euler
    IE, //Implicit Euler
    CN, //Crank-Nicholson
    PI  //Point-Implicit
  };

  /* Pure Implicit / Explicit / Defined Mixture */
  template<Integrator I>
  void integrate(double dt, Eigen::VectorXd& val, Eigen::SparseMatrix<double>& mat){
    std::cout<<"Please specify an integration method."<<std::endl;
  }

  /* Pure Implicit / Explicit / Defined Mixture + Source Term */
  template<Integrator I>
  void integrate(double dt, Eigen::VectorXd& val, Eigen::SparseMatrix<double>& mat, Eigen::VectorXd& source){
    std::cout<<"Please specify an integration method."<<std::endl;
  }

  /* Point Implicit / Arbitrary Mixture */
  template<Integrator I>
  void integrate(double dt, Eigen::VectorXd& val, Eigen::SparseMatrix<double>& mat_E, Eigen::SparseMatrix<double>& mat_I){
    std::cout<<"Please specify an integration method."<<std::endl;
  }

  template<Integrator I>
  void integrate(double dt, Eigen::VectorXd& val, Eigen::SparseMatrix<double>& mat_E, Eigen::SparseMatrix<double>& mat_I, Eigen::VectorXd& source){
    std::cout<<"Please specify an integration method."<<std::endl;
  }

  /* Explicit Euler Integrator - Fully Explicit */

  template<>
  void integrate<EE>(double dt, Eigen::VectorXd& val, Eigen::SparseMatrix<double>& mat){
    val.noalias() = (am::sparseIdentity()+dt*mat)*val;
  }

  template<>
  void integrate<EE>(double dt, Eigen::VectorXd& val, Eigen::SparseMatrix<double>& mat, Eigen::VectorXd& source){
    val.noalias() = (am::sparseIdentity()+dt*mat)*val + dt*source;
  }

  /* Implicit Euler Integrator - Fully Implicit */

  template<>
  void integrate<IE>(double dt, Eigen::VectorXd& val, Eigen::SparseMatrix<double>& mat){
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(am::sparseIdentity()-dt*mat);
    val.noalias() = solver.solve(val);
  }

  template<>
  void integrate<IE>(double dt, Eigen::VectorXd& val, Eigen::SparseMatrix<double>& mat, Eigen::VectorXd& source){
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(am::sparseIdentity()-dt*mat);
    val.noalias() = solver.solve(val + dt*source);
  }

  /* Crank-Nicholson Integrator - Semi Explicit / Implicit */

  template<>
  void integrate<CN>(double dt, Eigen::VectorXd& val, Eigen::SparseMatrix<double>& mat){
    integrate<EE>(0.5*dt, val, mat); //Perform Explicit Half-Step
    integrate<IE>(0.5*dt, val, mat); //Perform Implicit Half-Step
  }

  template<>
  void integrate<CN>(double dt, Eigen::VectorXd& val, Eigen::SparseMatrix<double>& mat, Eigen::VectorXd& source){
    integrate<EE>(0.5*dt, val, mat, source); //Perform Explicit Half-Step
    integrate<IE>(0.5*dt, val, mat, source); //Perform Implicit Half-Step
  }

  /* Point-Implicit Integrator */

  template<>
  void integrate<PI>(double dt, Eigen::VectorXd& val, Eigen::SparseMatrix<double>& mat_E, Eigen::SparseMatrix<double>& mat_I){
    integrate<EE>(dt, val, mat_E); //Perform Explicit Integration
    integrate<IE>(dt, val, mat_I); //Perform Implicit Half-Step
  }

  template<>
  void integrate<PI>(double dt, Eigen::VectorXd& val, Eigen::SparseMatrix<double>& mat_E, Eigen::SparseMatrix<double>& mat_I, Eigen::VectorXd& source){
    integrate<EE>(dt, val, mat_E, source); //Perform Explicit Integration
    integrate<IE>(dt, val, mat_I, source); //Perform Implicit Half-Step
  }

  /* Feeding 0.5 Mat for Both is equivalent to CN! */
  /* You don't even need to feed 0.5 in theory! */
};
