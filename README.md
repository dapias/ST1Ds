# Numerical tests for ergodicity

We provide the code that supports the numerical tests reported in the manuscript *Ergodicity of one--dimensional systems coupled to the logistic thermostat*.

It is organized as follows:

The `src` folder contains the files where the equations of motion are defined for a general potential together with the needed functions to calculate the Lyapunov spectrum and the Poincaré sections.

The `simulation` folder contains the main files `runsimulations.jl` and `poincareintegration.jl`. The first one generates a data file (with extension *.hdf5*) which contains the Lyapunov spectra for a certain number of random initial conditions specified in the file `parameterssimulation.yaml` (where the rest of relevant conditions can be specified). The second file integrates an arbitrary initial condition for a very long time, specified in `parametersintegration.yaml` and generates a data file which stores the integrated trajectory. Both data files are generated in a new folder that automatically will be createdcalled `data`.

Finally, the folder `analysis` contains the files needed to generate the kind of data reported in the article that may be imported and directly plotted without further treatment.

## Requirements

**General**
- Julia. 
**Packages**
- PyCall
- HDF5
- YAML
- ODE
- Cubature
- ForwardDiff

To add a package type the following command in the Julia REPL:

```
julia> Pkg.add("PackageName")
```

## Authors

Diego Tapias (Facultad de Ciencias, UNAM) diego.tapias@nucleares.unam.mx

David Sanders (Facultad de Ciencias, UNAM) dpsanders@ciencias.unam.mx



