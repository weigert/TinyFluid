class Input{
public:
  bool quit = false;
  bool paused = true;
  bool trigger = false;

  void handle();
  SDL_Event event;
};
