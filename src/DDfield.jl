using ForwardDiff
import ForwardDiff.derivative

type Thermostat{F<:Function}
    name::String
    Q::Float64
    distribution::F
end

type Potential{F<:Function}
    name::String
    f::F
end

function force(x::Float64, potential::Potential)
    f(y) = -derivative(potential.f,y)
    f(x)
end

function friction(z::Float64, thermo::Thermostat)
    g = y -> derivative(thermo.distribution,y)
    g(z)/thermo.distribution(z)
end

"""
Returns the derivative of the force given a potential,
i.e. the negative second derivative of the potential
"""
function forcederivative(potential::Potential)
    force(x) = -derivative(potential.f,x)
    fprime(x) = derivative(force,x)
end

    
"""
Equations of motion fot a general potential coupled to the logistic thermostat
"""
function DDfield(r::Vector{Float64}, potential::Potential, beta::Float64, thermo::Thermostat)

    (q, p, z) = r

    f = x-> force(x, potential)
    
    dq_dt = p
    dp_dt = f(q) + friction(z,thermo)*p/beta
    dz_dt = p^2. - 1.0/beta

    [dq_dt, dp_dt, dz_dt]

end




"""
    jacobian(r, potential, beta, Q)

Returns the jacobian of the DDfield evaluated at the point `r` for a given `potential` and some values of `beta` and `Q`
"""
function jacobian(r::Vector{Float64}, potential::Potential, beta::Float64, thermo::Thermostat)

    (q,p,z) = r
    
    fprime = forcederivative(potential)

    distprime(x) = derivative(thermo.distribution,x)
    distprime2(x) = derivative(distprime,x)
    
    
    J = [0. 1. 0.; fprime(q)  friction(z, thermo)/beta p/beta*(distprime2(z)/thermo.distribution(z) - friction(z, thermo)^2.); 0. 2.0*p 0.]

end
    
    
"""
    variationalDDfield(r, potential, beta, Q)

Evaluates the variational field associated to the DDfield at the extended vector `r` (q,p,z  and 9 more entries representing the coordinates of three tangent vectors) for a given `potential` and some values of `beta` and `Q`
"""
function variationalDDfield(r_and_phi::Vector{Float64}, potential::Potential, beta::Float64, thermo::Thermostat)

    v = DDfield(r_and_phi[1:3], potential, beta, thermo)
    J = jacobian(r_and_phi[1:3], potential, beta, thermo)
    rmatrix = reshape(r_and_phi[4:end],3,3)
    DPhi = J*rmatrix'

    return append!(v, DPhi'[:])
end    


function rungeK(x_init::Vector{Float64},field::Function,stepsize::Float64)
  state = x_init
  k1 = field(state)
  k2 = field(state + .5*stepsize*k1)
  k3 = field(state + .5*stepsize*k2)
  k4 = field(state + stepsize*k3)
  phi = (1/6.)*k1 + (1/3.)*k2 + (1/3.)*k3 + (1/6.)*k4
  state = state+stepsize*phi
  return state
end

"""
Solver based on the 4th order Runge-Kutta integrator
"""
function flowRK(field::Function, r0::Vector{Float64},dt::Float64, tfinal::Float64, potential::Potential, beta::Float64, thermo::Thermostat)

    t = 0.0:dt:tfinal
    pos = copy(r0)

    function extendedfield(r::Vector{Float64})
        field(r, potential, beta, thermo)
    end

    N = length(t) - 1
    for i in 1:N
        pos = rungeK(pos, extendedfield, dt)
    end

    return pos
end
