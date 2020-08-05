#include <iostream>

#define SIZE 50

#include "../../TinyFluid.h"
#include "flowtest.h"

int main(int argc, char* argv[]){

  /* Setup the Renderer */
  View view;
  if(!view.setup()){
    std::cout<<"Failed to setup view class."<<std::endl;
    return 0;
  }

  /* Initialize the Field Class */
  Field field;
  Input input;

  /* Initialize topography to some height-map */
  field.initialize();

  /* Render Loop */
  while(!input.quit){
    //Check for Quit
    while( SDL_PollEvent( &input.event ) != 0 ) input.handle();
    if(input.quit) break;

    if(!input.paused){
      field.timestep();
    }

    float _min = field.P.minCoeff();
    float _max = field.P.maxCoeff();

    //Drawing Logic
    view.render([&](){

      //Do this for all squares!
      for(int i = 0; i < SIZE*SIZE; i++){
        //Get the element Position
        glm::vec2 _pos = alg::pos(i);

        double P = (field.P(i) - _min)/(_max-_min);
        view.drawPixel(_pos, color::bezier(P, color::nebula), 1.0);

        //Velocity Field
        //double V = sqrt(pow(field.vX(i), 2) + pow(field.vY(i), 2))/2.5;
        //view.drawPixel(_pos, color::bezier(V, color::nebula), 1.0);

        //Draw the Quiver
        view.drawLine(_pos, glm::vec2(field.vX(i), field.vY(i)));
      }

		});
  }

  //Finished!
  return 0;
}
