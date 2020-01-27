class Input{
public:
  bool quit = false;
  bool paused = true;

  void handle();
  SDL_Event event;
};
