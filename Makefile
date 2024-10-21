CC = g++
CFLAGS = -std=c++11

SRCS = src/Ball.cpp src/Flow.cpp src/Simulation.cpp src/Collision.cpp src/main.cpp
OBJS = $(SRCS:.cpp=.o)

all: simulation

simulation: $(OBJS)
	$(CC) $(CFLAGS) -o simulation $(OBJS)

clean:
	rm -f *.o simulation
