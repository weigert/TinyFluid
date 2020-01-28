typedef Eigen::Triplet<double> triplet;
#include <math.h>

namespace am{

  //Convert Index to Coordinates and Back
  glm::ivec2 pos(int i){
    return glm::ivec2(i/SIZE, i%SIZE);
  }

  int index(glm::vec2 pos){
    return pos.x*SIZE + pos.y;
  }

  int index(int x, int y){
    return x*SIZE + y;
  }

  //Modulus Operator
  glm::vec2 operator%(const glm::vec2& rhs, const int x){
    return glm::mod(rhs, glm::vec2(x));
  }

  glm::ivec2 operator%(const glm::ivec2& rhs, const int x){
    return glm::mod(glm::vec2(rhs.x, rhs.y), glm::vec2(x));
  }

  std::ostream& operator<<(std::ostream& out, const glm::vec2 vec){
    std::cout<<"X: "<<vec.x<<" Y: "<<vec.y;
    return out;
  }


  Eigen::SparseMatrix<double> sparseIdentity(){
    Eigen::SparseMatrix<double> I(SIZE*SIZE, SIZE*SIZE);
    std::vector<triplet> list;
    for(int i = 0; i < SIZE*SIZE; i++){
      list.push_back(triplet(i, i, 1.0));
    }
    I.setFromTriplets(list.begin(), list.end());
    return I;
  }

  Eigen::SparseMatrix<double> sparseDiagonalize(Eigen::VectorXd &vec){
    Eigen::SparseMatrix<double> I(SIZE*SIZE, SIZE*SIZE);
    std::vector<triplet> list;
    for(int i = 0; i < SIZE*SIZE; i++){
      list.push_back(triplet(i, i, vec(i)));
    }
    I.setFromTriplets(list.begin(), list.end());
    return I;
  }

  Eigen::VectorXd flatGaussian(glm::vec2 mean, double stddev){
    Eigen::VectorXd gauss = Eigen::ArrayXd::Zero(SIZE*SIZE);
    for(int i = 0; i < SIZE*SIZE; i++){
      glm::vec2 _pos = pos(i);
      gauss(i) = exp(-(pow(_pos.x - mean.x, 2) + pow(_pos.y - mean.y, 2))/stddev);
    }
    return gauss;
  }


  /*
    These Operators could theoretically be rewritten for arbitrary order.
    They should be rewritten in general form for arbitrary dimension (not too hard)
  */

  /*
    Two Types of Operators:
      -> Interpolators
      -> Finite Differences

    Then all subsequent operators are simply combinations of the two.
  */

  /*
    Lagrange Interpolator:

    What is the value at a specific point, given that we know the surrounding points?

    I know two surround points -> Compute the linear interpolation polynomial

    Then express the value at our specific point as the weighted sum.

    These weights are the entry in our matrix.
  */

  std::vector<double> lagrange(std::initializer_list<int> shift){
    std::vector<double> weights;

    /* This is equivalent to a matrix inversion problem. */

    return weights;
  }

  /*

    Finite Difference Approximator:

    Define an arbitrary order of derivation, and an arbitrary number of point offsets.

    Then return the matrix that will give that approximation.

    Construct the Taylor Matrix and Invert it (or use some formula!!)

  */

  Eigen::SparseMatrix<double> FD(std::initializer_list<int> points, int order){
    Eigen::SparseMatrix<double> M(SIZE*SIZE, SIZE*SIZE);
    return M;
  }

  //All other Matrices are Derived from the FD Method and the Lagrange Interpolator.
  //Then, you can add a lagrange surface interpolator.
  //Then, you can add multidimensional Finite Differences.


  //Derived Functions

  //Mid-Point Calculation Matrix (For Surface Integrals)
  Eigen::SparseMatrix<double> getTrapezeX(int back){
    //Initialize Sparse Matrix
    Eigen::SparseMatrix<double> M(SIZE*SIZE, SIZE*SIZE);
    std::vector<triplet> list;

    //Loop over all Indices
    for(int i = 0; i < SIZE*SIZE; i++){
      //Position of this Element
      glm::ivec2 _pos = pos(i);
      //glm::ivec2 _newposA = (_pos + glm::ivec2(back, 1) + SIZE)%SIZE;
      //glm::ivec2 _newposB = (_pos + glm::ivec2(back,-1) + SIZE)%SIZE;
      glm::ivec2 _newpos = (_pos + glm::ivec2(back, 0) + SIZE)%SIZE;

      //Add Triplets
      list.push_back(triplet(i, i, 0.5));
      list.push_back(triplet(i, index(_newpos), 0.5));
      //list.push_back(triplet(i, index(_newposA), 0.25));
      //list.push_back(triplet(i, index(_newposB), 0.25));
    }

    //Construct and Return M
    M.setFromTriplets(list.begin(), list.end());
    return M;
  }

  Eigen::SparseMatrix<double> getTrapezeY(int back){
    //Initialize Sparse Matrix
    Eigen::SparseMatrix<double> M(SIZE*SIZE, SIZE*SIZE);
    std::vector<triplet> list;

    //Loop over all Indices
    for(int i = 0; i < SIZE*SIZE; i++){
      //Position of this Element
      glm::ivec2 _pos = pos(i);
      //glm::ivec2 _newposA = (_pos + glm::ivec2(SIZE+1, SIZE+back))%SIZE;
      //glm::ivec2 _newposB = (_pos + glm::ivec2(SIZE-1, SIZE+back))%SIZE;
      glm::ivec2 _newpos = (_pos + glm::ivec2(0, back) + SIZE)%SIZE;

      //Add Triplets
      //list.push_back(triplet(i, i, 0.5));
      //list.push_back(triplet(i, index(_newposA), 0.25));
      //list.push_back(triplet(i, index(_newposB), 0.25));
      list.push_back(triplet(i, i, 0.5));
      list.push_back(triplet(i, index(_newpos), 0.5));
    }

    //Construct and Return M
    M.setFromTriplets(list.begin(), list.end());
    return M;
  }

  Eigen::SparseMatrix<double> getDiffX(int back){
    //Initialize Sparse Matrix
    Eigen::SparseMatrix<double> M(SIZE*SIZE, SIZE*SIZE);
    std::vector<triplet> list;

    for(int i = 0; i < SIZE*SIZE; i++){
      //Position of this Element
      glm::ivec2 _pos = pos(i);
      glm::ivec2 _newpos = (_pos + glm::ivec2(SIZE+back, 0))%SIZE;

      //Add Triplets
      list.push_back(triplet(i, i, -1.0));
      list.push_back(triplet(i, index(_newpos), 1.0));
    }

    //Construct and Return M
    M.setFromTriplets(list.begin(), list.end());
    return M;
  }

  /* High Order Approximation */
  Eigen::SparseMatrix<double> getDiffY(int back){
    //Initialize Sparse Matrix
    Eigen::SparseMatrix<double> M(SIZE*SIZE, SIZE*SIZE);
    std::vector<triplet> list;

    for(int i = 0; i < SIZE*SIZE; i++){
      //Position of this Element
      glm::ivec2 _pos = pos(i);
      glm::ivec2 _newpos = (_pos + glm::ivec2(0, SIZE+back))%SIZE;

      //Add Triplets
      list.push_back(triplet(i, i, -1.0));
      list.push_back(triplet(i, index(_newpos), 1.0));
    }
    //Construct and Return M
    M.setFromTriplets(list.begin(), list.end());
    return M;
  }

  /* High Order Approximation */
  Eigen::SparseMatrix<double> gradX(){
    //Initialize Sparse Matrix
    Eigen::SparseMatrix<double> M(SIZE*SIZE, SIZE*SIZE);
    std::vector<triplet> list;

    for(int i = 0; i < SIZE*SIZE; i++){
      //Position of this Element
      glm::ivec2 _pos = pos(i);
      glm::ivec2 _newposB2 = (_pos + glm::ivec2(-2, 0) + SIZE)%SIZE;
      glm::ivec2 _newposB1 = (_pos + glm::ivec2(-1, 0) + SIZE)%SIZE;
      glm::ivec2 _newposF1 = (_pos + glm::ivec2( 1, 0) + SIZE)%SIZE;
      glm::ivec2 _newposF2 = (_pos + glm::ivec2( 2, 0) + SIZE)%SIZE;

      //Add Triplets
      list.push_back(triplet(i, index(_newposB2),  1.0/12.0));
      list.push_back(triplet(i, index(_newposB1), -8.0/12.0));
      list.push_back(triplet(i, index(_newposF1),  8.0/12.0));
      list.push_back(triplet(i, index(_newposF2), -1.0/12.0));
    }
    //Construct and Return M
    M.setFromTriplets(list.begin(), list.end());
    return M;
  }

  /* High Order Approximation */
  Eigen::SparseMatrix<double> gradY(){
    //Initialize Sparse Matrix
    Eigen::SparseMatrix<double> M(SIZE*SIZE, SIZE*SIZE);
    std::vector<triplet> list;

    for(int i = 0; i < SIZE*SIZE; i++){
      //Position of this Element
      glm::ivec2 _pos = pos(i);
      glm::ivec2 _newposB2 = (_pos + glm::ivec2( 0,-2) + SIZE)%SIZE;
      glm::ivec2 _newposB1 = (_pos + glm::ivec2( 0,-1) + SIZE)%SIZE;
      glm::ivec2 _newposF1 = (_pos + glm::ivec2( 0, 1) + SIZE)%SIZE;
      glm::ivec2 _newposF2 = (_pos + glm::ivec2( 0, 2) + SIZE)%SIZE;

      //Add Triplets
      list.push_back(triplet(i, index(_newposB2),  1.0/12.0));
      list.push_back(triplet(i, index(_newposB1), -8.0/12.0));
      list.push_back(triplet(i, index(_newposF1),  8.0/12.0));
      list.push_back(triplet(i, index(_newposF2), -1.0/12.0));
    }
    //Construct and Return M
    M.setFromTriplets(list.begin(), list.end());
    return M;
  }

  /* High Order Approximation */
  Eigen::SparseMatrix<double> laplaceX(){
    //Initialize Sparse Matrix
    Eigen::SparseMatrix<double> M(SIZE*SIZE, SIZE*SIZE);
    std::vector<triplet> list;

    for(int i = 0; i < SIZE*SIZE; i++){
      //Position of this Element
      glm::ivec2 _pos = pos(i);
      glm::ivec2 _newposB2 = (_pos + glm::ivec2(-2, 0) + SIZE)%SIZE;
      glm::ivec2 _newposB1 = (_pos + glm::ivec2(-1, 0) + SIZE)%SIZE;
      glm::ivec2 _newposF1 = (_pos + glm::ivec2( 1, 0) + SIZE)%SIZE;
      glm::ivec2 _newposF2 = (_pos + glm::ivec2( 2, 0) + SIZE)%SIZE;

      //Add Triplets
      list.push_back(triplet(i, index(_newposB2), -1.0/12.0));
      list.push_back(triplet(i, index(_newposB1), 16.0/12.0));
      list.push_back(triplet(i, i,  -30.0/12.0));
      list.push_back(triplet(i, index(_newposF1), 16.0/12.0));
      list.push_back(triplet(i, index(_newposF2), -1.0/12.0));
    }
    //Construct and Return M
    M.setFromTriplets(list.begin(), list.end());
    return M;
  }

  /* High Order Approximation */
  Eigen::SparseMatrix<double> laplaceY(){
    //Initialize Sparse Matrix
    Eigen::SparseMatrix<double> M(SIZE*SIZE, SIZE*SIZE);
    std::vector<triplet> list;

    for(int i = 0; i < SIZE*SIZE; i++){
      //Position of this Element
      glm::ivec2 _pos = pos(i);
      glm::ivec2 _newposB2 = (_pos + glm::ivec2( 0,-2) + SIZE)%SIZE;
      glm::ivec2 _newposB1 = (_pos + glm::ivec2( 0,-1) + SIZE)%SIZE;
      glm::ivec2 _newposF1 = (_pos + glm::ivec2( 0, 1) + SIZE)%SIZE;
      glm::ivec2 _newposF2 = (_pos + glm::ivec2( 0, 2) + SIZE)%SIZE;

      //Add Triplets
      list.push_back(triplet(i, index(_newposB2), -1.0/12.0));
      list.push_back(triplet(i, index(_newposB1), 16.0/12.0));
      list.push_back(triplet(i, i,  -30.0/12.0));
      list.push_back(triplet(i, index(_newposF1), 16.0/12.0));
      list.push_back(triplet(i, index(_newposF2), -1.0/12.0));
    }
    //Construct and Return M
    M.setFromTriplets(list.begin(), list.end());
    return M;
  }

}; //End of Namespace
