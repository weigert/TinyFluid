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
	int ratiox = SCREEN_WIDTH / SIZE;
	int ratioy = SCREEN_HEIGHT / SIZE;
  SDL_Rect rect{ratiox*pos.x, ratioy*pos.y, ratiox, ratioy};
  SDL_SetRenderDrawColor(gRenderer, color.x, color.y, color.z, 255);
  SDL_RenderFillRect(gRenderer, &rect);
}

void View::drawLine(float x, float y, float dx, float dy){
	int ratiox = SCREEN_WIDTH / SIZE;
	int ratioy = SCREEN_HEIGHT / SIZE;

	/* I need Direction AND Intensity */
	SDL_SetRenderDrawColor(gRenderer, 0, 0, 0, 255);
	SDL_RenderDrawLine(gRenderer, ratiox*(x+0.5-dx), ratioy*(y+0.5-dy), ratiox*(x+0.5+dx), ratioy*(y+0.5+dy));
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
