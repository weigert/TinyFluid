/*
  Physical Source Terms

  Desirable Terms:
  - Henry's Law for Vapor equilibrium concentration
  - Raoults Law for Determining Phase Fraction!
  - Thin Film Mass Transfer Law
  - Sun Heating (Flat or Radial)
  -
*/

namespace source{

  Eigen::VectorXd Psat(Eigen::VectorXd& P, Eigen::VectorXd& T){
    Eigen::VectorXd _psat = Eigen::ArrayXd::Ones(SIZE*SIZE);
    for(int i = 0; i < SIZE*SIZE; i++)
      _psat(i) = pow(10.0, 8.07131 - 1730.63/(233.426+T(i))) * 101325.0/760.0;
    return _psat;
  }


  //How do clouds look?
  Eigen::VectorXd CLOUD(Eigen::VectorXd& H, Eigen::VectorXd& P, Eigen::VectorXd& T, double thresh){
    Eigen::VectorXd E = Eigen::ArrayXd::Ones(SIZE*SIZE);
    Eigen::VectorXd MTM = thresh*E-(Psat(P, T)-H.cwiseProduct(P)).cwiseMin(thresh);  //Too large a discrepancy!
    return MTM/thresh*0.9; //It is possible to be "oversaturated"
  }

  Eigen::VectorXd TSOURCE(Eigen::VectorXd &B, Eigen::VectorXd& height, double s, Eigen::VectorXd& H, Eigen::VectorXd& P, Eigen::VectorXd& T, Eigen::SparseMatrix<double> &vX_MAT, Eigen::SparseMatrix<double> &vY_MAT){

    Eigen::VectorXd E = Eigen::ArrayXd::Ones(SIZE*SIZE);

    //Evaporation Cooling (MORE!)
    double m = 0.001;
    Eigen::VectorXd MTD = (Psat(P, T)-H.cwiseProduct(P)).cwiseMax(0.0);
    Eigen::VectorXd Z = -50000.0*m*MTD.cwiseProduct(s*E-height.cwiseMin(s)).cwiseQuotient(Psat(P, T));

    //Rain Removal
    double thresh = 75.0; //This is the oversaturation capacity
    double k = 0.0001;
    Eigen::VectorXd MTM = -E*thresh+(Psat(P, T)-H.cwiseProduct(P)).cwiseMin(thresh);  //Too large a discrepancy!
    Z += 1000.0*k*MTM;

    //Boundary Volume Not Height Difference... (Increase proportionally to the temperature...)
    Z -= 50.0*(0.5*vX_MAT*(PDE::GXF+PDE::GXB)*height.cwiseMax(s) + 0.5*vY_MAT*(PDE::GYF+PDE::GYB)*height.cwiseMax(s));

    //Radiation Transfer (Get Albedo and Radiation Balance Steady State T)
    Eigen::VectorXd C = CLOUD(H, P, T, 250.0);
    Z += 0.1*(10.0*E + 15*E.cwiseProduct(shape::flatGaussian(glm::vec2(SIZE/2.0, SIZE/2.0), glm::vec2(SIZE*15.0, SIZE*2.0)))-T).cwiseProduct(E-C);

    return Z;
  }

  Eigen::VectorXd HSOURCE(Eigen::VectorXd& B, Eigen::VectorXd& height, double s, Eigen::VectorXd& H, Eigen::VectorXd& P, Eigen::VectorXd& T){

    Eigen::VectorXd E = Eigen::ArrayXd::Ones(SIZE*SIZE);

    //Mass Transfer Addition..
    //Evaporation Mass Transfer (k = approach rate to th equilibrium concentration)
    double m = 0.001;
    Eigen::VectorXd MTD = (Psat(P, T)-H.cwiseProduct(P)).cwiseMax(0.1);
    Eigen::VectorXd Z = m*MTD.cwiseProduct(B).cwiseQuotient(Psat(P, T));

    //Rain Removal
    double thresh = 75.0; //This is the oversaturation capacity!
    double k = 0.0001;
    Eigen::VectorXd MTM = -E*thresh+(Psat(P, T)-H.cwiseProduct(P)).cwiseMin(thresh);  //Too large a discrepancy!
    Z += k*MTM;

    return Z;
  }
};
