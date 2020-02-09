/*
Author: Nicholas McDonald

A bunch of helper functions that can arrays that represent a shape in N dimensions.

Example:
  Sphere / Circle as a boolean array in 3D / 2D ...
  Gaussians in 1D / 2D / 3D Centered around a point, filled with floats...
  etc.

This is useful because arrays need to be flattened / unflattened.
*/


namespace shape{
  //First Shape Function!
  Eigen::VectorXd flatGaussian(glm::vec2 mean, double stddev){
    Eigen::VectorXd gauss = Eigen::ArrayXd::Zero(SIZE*SIZE);
    for(int i = 0; i < SIZE*SIZE; i++){
      glm::vec2 _pos = alg::pos(i);
      gauss(i) = exp(-(pow(_pos.x - mean.x, 2) + pow(_pos.y - mean.y, 2))/stddev);
    }
    return gauss;
  }

  Eigen::VectorXd circle(glm::vec2 center, double rad){
    Eigen::VectorXd circ = Eigen::ArrayXd::Ones(SIZE*SIZE);
    for(int i = 0; i < SIZE*SIZE; i++){
      glm::vec2 _pos = alg::pos(i);
      if(sqrt(pow(_pos.x-center.x, 2) + pow(_pos.y-center.y, 2)) < rad)
        circ(i) = 0.0;
    }
    return circ;
  }
};
