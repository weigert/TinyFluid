################################################################################
# Climate Simulator
################################################################################
# Parameters
CC = g++ -std=c++14
FLAGS = -march=native -Wall -g -rdynamic -O3
LIBS	= -I/usr/local/include -L/usr/local/lib -lnoise -lSDL2 -pthread

# Build Targets
A_SRC = main.cpp
A_TARGET  = flowtest

# Make Rules
.PHONY: all
all: climate

.PHONY: climate
climate: $(A_SRC)
	$(CC) $(A_SRC) $(FLAGS) $(LIBS) -o $(A_TARGET)
