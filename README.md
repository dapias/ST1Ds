# Numerical tests for ergodicity of Singly 1D Thermostatted Systems (ST1Ds) 

We provide the code that supports the numerical tests reported in the manuscript [Ergodicity of one--dimensional systems coupled to the logistic thermostat]().

It is organized as follows:

The `src` folder contains the files where the equations of motion are defined for a general potential together with the needed functions to calculate the Lyapunov spectrum and the Poincaré sections.

The `simulation` folder contains the main file `runsimulations.jl`. This one generates a set of data files (with extension *.hdf5*) depending on the results asked. It may generate a file with the Lyapunov spectra for a certain number of random initial conditions specified in the file `parameters.jl` (where the rest of relevant conditions can be specified). It may also generate different files where the trajectory for a certain initial condition is saved. Both sets of data files are generated in a new folder that automatically will be created.

Finally, the folder `analysis` contains the files needed to generate the kind of results reported in the article that may be imported and directly plotted without further treatment. 

## Authors

Diego Tapias (Facultad de Ciencias, UNAM) diego.tapias@nucleares.unam.mx

David Sanders (Facultad de Ciencias, UNAM) dpsanders@ciencias.unam.mx

*2016*


