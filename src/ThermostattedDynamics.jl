__precompile__(true)

module ThermostattedDynamics

using ODE
using HDF5
using PyCall
using Cubature
using ForwardDiff
using StatsBase

import Distributions.Normal
import ForwardDiff.derivative

export Thermostat, Potential, Integrator, Parameters, lyapunov_exponents, trajectory, marginaldistributions, hellingerdistance, sections,
save_lyapunov, save_trajectory, section

const sm = PyNULL()

function __init__()
    copy!(sm, pyimport_conda("statsmodels.api", "statsmodels"))
end

include("DensityDynamics/randominitialcondition.jl")
include("DensityDynamics/DDtypes.jl")
include("DensityDynamics/DDmethods.jl")
include("DensityDynamics/lyapunov.jl")
include("DensityDynamics/trajectory.jl")
include("DensityDynamics/runsimulations.jl")
include("analysis/marginaldistributions.jl")
include("analysis/hellingerdistance.jl")
include("analysis/sections.jl")



end

