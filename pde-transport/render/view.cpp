#include "view.h"

/* Setup and Cleanup */
bool View::setup(){
  //Initialize SDL
	if( SDL_Init( SDL_INIT_VIDEO ) < 0 ) return false;

  //Create window
  gWindow = SDL_CreateWindow( WINDOW_NAME.c_str(), SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN );
  if( gWindow == NULL ) return false;

  //Prepare the Renderer
  gRenderer = SDL_CreateRenderer(gWindow, -1, SDL_RENDERER_ACCELERATED);
  SDL_SetRenderDrawBlendMode(gRenderer, SDL_BLENDMODE_BLEND);

  return true;
}

void View::cleanup(){
  //Destroy window
	SDL_DestroyWindow( gWindow );
	SDL_DestroyRenderer (gRenderer);

	//Quit SDL subsystems
	SDL_Quit();
}

/* Drawing Helpers */
void View::drawPixel(glm::ivec2 pos, glm::ivec3 color){
	/* Construct a Rect and Fill with Color at Position */
  SDL_Rect rect{10*pos.x, 10*pos.y, 10, 10};
  SDL_SetRenderDrawColor(gRenderer, color.x, color.y, color.z, 255);
  SDL_RenderFillRect(gRenderer, &rect);
}

void View::drawLine(float x, float y, float dx, float dy){
	float scale = 2.5;
	/* I need Direction AND Intensity */
	SDL_SetRenderDrawColor(gRenderer, 0, 0, 0, 255);
	SDL_RenderDrawLine(gRenderer, 10*x-dx*scale, 10*y-dy*scale, 10*x+dx*scale, 10*y+dy*scale);
}

template<typename F, typename... Args>
void View::render(F function, Args&&... args){
	//Clear the Window
	SDL_SetRenderDrawColor(gRenderer, 0, 0, 0, 255);
  SDL_RenderClear(gRenderer);

	//Call whatever the user has specified...
	function(args...);

	//Present the Information
	SDL_RenderPresent(gRenderer);
}
