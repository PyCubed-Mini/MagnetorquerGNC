## Simulator 

Contains code necessary to generate a simple satellite orbit around Earth. 
Updates both translational (position, velocity) states and rotational (attitude, angular velocity) states. 

Module exports several helpful functions, including:
  - rk4: Integrator function that takes in a current state and provides the next state 
  - IGRF13: Function that takes in current position and time and determines the magnetic field vector in inertial frame 
  - initialize_orbit: Randomly generates a valid orbit 
  - initialize_orbit_oe: Randomly generates a valid orbit, allowing for orbital elements to be specified.


There are also some (very simple) tests in the ```tests/``` folder that give an overview on how the functions can be called. 