include("./DDfield.jl")
include("./RungeKutta.jl")
include("./randominitialcondition.jl")

using ODE

function flowode45(field::Function, r0::Vector{Float64},dt::Float64, tfinal::Float64, potential::Function, beta::Float64, Q::Float64)

    t = 0.0:dt:tfinal

    function extendedfield(time, r::Vector{Float64})
        field(r, potential, beta, Q)
    end

    (t, pos) = ode45(extendedfield, r0, t, points=:specified)

end


