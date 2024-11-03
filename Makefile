Hello: Hello.cpp
	g++ Hello.cpp -o Hello

Balls: Main.cpp Balls.cpp Collision.cpp Configuration.cpp Flow.cpp Simulation.cpp Balls.h Collision.h Configure.h Flow.h Simulation.h 
	g++ Balls-in-poiseuille-fluids-with-oneself-collisions.cpp -o Balls
