Two types of simulation:

- Lyapunov spectrum for three potentials (Mexican Hat, Quartic oscillator and Harmonic oscillator) with different random 
initial conditions (**simulation 1**).

- Poincaré section for a particular long trajectory corresponding to one potential (**simulation 2**).

### Simulation 1

Run the file *simulations/runsimulations.jl*. Set the desired parameters in the file *simulations/parameterssimulation.yaml*. 
That file generates a HDF5 file that is saved in the folder *data/* with the lyapunov exponents and the initial conditions 
corresponding to the number of simulations set. Then you can proceed to group in a single data structure the results of a whole
group of simulations, to do that execute the file *analysis/lyapunovexponents.jl*. Finally you can group together the 
results of different groups of simulations by running the file *analysis/analysislyapunov.jl*. This script also 
generates and saves three histograms corresponding to the distribution of frequencies of each exponent.

### Simulation 2

Run the file *simulations/poincareintegration.jl*. Set the desired parameters in the file *simulations/parametersintegration.yaml*. 
That file generates a HDF5 file that is saved in the folder *poincaredata/* with the position dat $(p,q,z)$ corresponding to a 
particular trajectory integrated by the Dormand-Prince method for the parameters set. Then you can proceed to execute the file 
*simulations/plotpoincare.jl* which displays and saves the poincare section of the called file.



