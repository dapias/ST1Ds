module ThermostattedDynamics

using ODE
using HDF5
using PyCall
using Cubature
using ForwardDiff
using StatsBase

import Distributions.Normal

include("DensityDynamics/DDtypes.jl")
include("DensityDynamics/DDmethods.jl")
include("DensityDynamics/runsimulations.jl")
include("analysis/marginaldistributions.jl")
include("analysis/hellingerdistance.jl")
include("analysis/sections.jl")

export Thermostat, Potential, Integrator, Parameters, runsimulation, marginaldistributions, hellingerdistance, sections

end
