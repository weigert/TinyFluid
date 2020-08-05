#include <iostream>

#define SIZE 100

#include "../../TinyFluid.h"
#include "boundary.h"

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

  //Max and Min-Values! (Initial)
  //Pressure Field (Scaled)


  /* Render Loop */
  while(!input.quit){
    //Check for Quit
    while( SDL_PollEvent( &input.event ) != 0 ) input.handle();
    if(input.quit) break;

    if(!input.paused){
      timer::benchmark<std::chrono::microseconds>([&](){
        field.timestep();
      });
    }

    //Drawing Logic
    std::cout<<"Render ";
    timer::benchmark<std::chrono::microseconds>([&](){

    view.render([&](){

      //Do this for all squares!
      for(int i = 0; i < SIZE*SIZE; i++){
        //Get the element Position
        glm::vec2 _pos = alg::pos(i);

        double _V = sqrt(pow(field.vX(i), 2) + pow(field.vY(i), 2))/2.5;
        view.drawPixel(_pos, color::bezier(_V, color::nebula), 1.0);

        double _B = field.B(i);
        view.drawPixel(_pos, color::bezier(_B, color::bw), 1.0-_B);

        //double min = field.P.minCoeff();
        //double max = field.P.maxCoeff();
        //double P = (field.P(i) - min)/(max-min);
        //view.drawPixel(_pos, color::bezier(P, color::nebula), B);

        //Draw the Quiver
        view.drawLine(_pos, glm::vec2(field.vX(i), field.vY(i)));
      }
		});

    });
  }

  //Finished!
  return 0;
}
