#include <iostream>
#include <SDL2/SDL.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseLU>
#include <glm/glm.hpp>
#include <unistd.h>

#define SIZE 50

#include "../../pde-transport/pde-transport.h"
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

    //Drawing Logic
    view.render([&](){

      for(int i = 0; i < SIZE*SIZE; i++){
        //Get the element Position
        glm::vec2 _pos = am::pos(i);

        //Draw the Color
        //double k = abs(field.vX(i))+abs(field.vY(i));
        double k = abs(field.P(i));
        view.drawPixel(_pos, glm::ivec3(0.001*abs(k)));

        //Draw the Quiver
        view.drawLine(_pos.x, _pos.y, field.vX(i), field.vY(i));
      }

		});
  }

  //Finished!
  return 0;
}
