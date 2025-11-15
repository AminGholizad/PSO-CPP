# Particle swarm optimization (PSO)
Particle swarm optimization algorithm ([PSO](https://en.wikipedia.org/wiki/Particle_swarm_optimization)) for a minimization problem. In this project, nonlinear constraints are implemented as infeasible solutions.
This project is using C++23 features.
Constraints and objective in one function.

# Features

1. Template class for particles with number of variables as the template parameter is used to initialize the arrays of the needed size.
2. Template function for pso with number of variables and swarm size is as the template parameters is used to initialize arrays accordingly
3. The cost function signiture should be one array/vector of doubles as input and the output is a cost type of an array of doubles representing objective and a double representing infeasablity of the problem.
4. The Mutation is used to avoid local minima.

# todo
[] More explanation of the functionality
