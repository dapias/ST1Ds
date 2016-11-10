include("./DDfield.jl")
include("./RungeKutta.jl")
include("./randominitialcondition.jl")

using ODE

"""
    flowode45(densityfield, r0, dt, tfinal, potential, beta, Q)

Function that integrates a trajectory for a given initial condition `r0`, given a vector DD field associated to  a certain potential, using the Dormand-Pince 45 integrator
"""
function flowode45(densityfield::Function, r0::Vector{Float64},dt::Float64, tfinal::Float64, potential::Function, beta::Float64, Q::Float64)

    t = 0.0:dt:tfinal

    function extendedfield(time, r::Vector{Float64})
        densityfield(r, potential, beta, Q)
    end

    (t, pos) = ode45(extendedfield, r0, t, points=:specified)

end


