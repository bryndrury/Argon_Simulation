# Argon_Simulation

This is a personal project, use at your own risk.

The main computation for this project is done through c++ and the code should work standalone.
There is a python script called sim_.py that is meant to help automate running multiple simulations at once.
Just populate the lists N and L in the python file sim_.py with the number of particles and the size of the container repectivly and run.

The conditions.py creates random initial conditions for the particles.
The animate.py uses matplotlib to create a gif of the particles motion.
The animation uses a lot of memory so you might need to lower the dpi (at th bottom of the file) or increase the Tsave value in the sim_.py file. 

**Warning** The program runs all the simulations and animations at once so limit the size of the N and L lists if your computer is but a mere mortal.
