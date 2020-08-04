#include <iostream>

#define SIZE 100

#include "../../procfluid/procfluid.h"
#include "terrain.h"

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
//  srand(time(NULL));
  field.SEED = 100;//rand()%10000;
  field.initialize();

  /* Render Loop */
  bool exit = false;
  while(!input.quit && !exit){


    //Check for Quit
    while( SDL_PollEvent( &input.event ) != 0 ) input.handle();
    if(input.quit) exit = true;

    timer::benchmark<std::chrono::microseconds>([&](){

    if(!input.paused){
        field.timestep();
    }

    });

    if(input.trigger){
      field.g = glm::vec2(0.0, 0.0);
      input.trigger = false;
      std::cout<<"Triggered"<<std::endl;
    }

    Eigen::VectorXd clouds;
    Eigen::VectorXd clouds2;

    //Drawing Logic
    view.render([&](){

      clouds = source::CLOUD(field.humidity, field.P, field.temperature, 100.0);
      clouds2 = source::CLOUD(field.humidity, field.P, field.temperature, 150.0);

      //Do this for all squares!
      for(int i = 0; i < SIZE*SIZE; i++){
        //Get the element Position
        glm::vec2 _pos = alg::pos(i);

        double H = field.height(i);
        if(H < field.sealevel) view.drawPixel(_pos, color::bezier(H/field.sealevel*0.5, color::ocean), 1.0);
        else view.drawPixel(_pos, color::bezier((H-field.sealevel)/(1.0-field.sealevel), color::land), 1.0);

        if(input.screen == 1){
          double min = field.temperature.minCoeff();
          double max = field.temperature.maxCoeff();
          std::cout<<"MIN: "<<min<<" MAX: "<<max<<std::endl;
          double T = (field.temperature(i)-min)/(max-min);
          view.drawPixel(_pos, color::red, T);
        }
        else if(input.screen == 2){
          double H = field.humidity(i);
          view.drawPixel(_pos, color::white, 50*H);
        }
        else if(input.screen == 3){
          double C = (clouds(i)>1.0)?1.0:clouds(i);
          double C2 = (clouds2(i)>1.0)?1.0:clouds2(i);
          view.drawPixel(_pos, color::white, 0.8*C);
          view.drawPixel(_pos, color::white, 0.1*C2);
        }
        else if(input.screen == 4){
          double V = sqrt(pow(field.vX(i), 2) + pow(field.vY(i), 2));
          view.drawPixel(_pos, color::bezier(V, color::nebula), 0.8);
        }
        else if(input.screen == 5){
          double min = field.P.minCoeff();
          double max = field.P.maxCoeff();
          double P = (field.P(i) - min)/(max-min);
          view.drawPixel(_pos, color::bezier(P, color::bw), 0.6);
        }

        //Draw the Quiver
        view.drawLine(_pos, glm::vec2(field.vX(i), field.vY(i)));
      }
		});

  }

  //Finished!
  return 0;
}
