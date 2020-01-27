#define SCREEN_WIDTH 500
#define SCREEN_HEIGHT 500

#include <string>
#include <SDL2/SDL.h>

//Master Class
class View {
 public:
   std::string WINDOW_NAME = "Climate Simulator";

   //Rendering Variables
   SDL_Window* gWindow;
   SDL_Renderer* gRenderer;

   //Setup and Cleanup
   bool setup();
   void cleanup();

   //Drawing Helpers
   void drawPixel(glm::ivec2 pos, glm::ivec3 color);
   void drawLine(float x, float y, float dx, float dy);

   //Renderer
   template<typename F, typename... Args>
   void render(F function, Args&&... args);
};
