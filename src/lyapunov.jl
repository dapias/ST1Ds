include("./ThermostatModel.jl")
include("./RungeKutta.jl")
include("./initialconditions.jl")

using ThermostatModel
using SymPy
using ForwardDiff
import ForwardDiff.derivative


function thermostat(thermotype::String, Q::Float64, T::Float64)
    thermomodel = eval(parse(thermotype))
    thermo = thermomodel(Q, 1/T)
end    

function DDfield(r::Vector{Float64}, potential::Function, thermo::Thermostat)

    (q, p, z) = r

    force(x::Float64) = -derivative(potential,x)

    
    dq_dt = p
    dp_dt = force(q) + p*friction(thermo, z)/thermo.beta
    dz_dt = p^2. - 1.0/thermo.beta

    [dq_dt, dp_dt, dz_dt]

end

function jacobian(r_and_phi::Vector{Float64}, potential::Function, thermo::Thermostat)

    x = Sym("x")
    fprime = -diff(diff(potential(x)))
    
    J = [0. 1. 0.; fprime(x=>r_and_phi[1]) friction(thermo, r_and_phi[3])/thermo.beta difffriction(thermo, r_and_phi[3])/thermo.beta; 0. 2.0*r_and_phi[2] 0.]

end
    
    

function variationalDDfield(r_and_phi::Vector{Float64}, potential::Function, thermo::Thermostat)

    r = DDfield(r_and_phi[1:3], potential, thermo)

    J = jacobian(r_and_phi, potential, thermo)
   
    rmatrix = reshape(r_and_phi[4:end],3,3)
    DPhi = J*rmatrix'
    

    return append!(r, DPhi'[:])

end    

function flow(field::Function, r0::Vector{Float64},dt::Float64, tfinal::Float64, potential::Function, thermo::Thermostat)

    t = 0.0:dt:tfinal
    pos = copy(r0)

    function extendedfield(r::Vector{Float64})
        field(r, potential, thermo)
    end

    N = length(t) - 1
    for i in 1:N
        pos = rungeK(pos, extendedfield, dt)
    end

    return pos
end
    
function lyapunovspectra(field::Function, r::Vector{Float64}, dt::Float64, dtsampling::Float64, nsteps::Int64, potential::Function, thermo::Thermostat)
    w = eye(3)
    norm1 = Float64[]
    norm2 = Float64[]
    norm3 = Float64[]
    for i in 1:nsteps
        pos = flow(variationalDDfield, r, dt, dtsampling, potential, thermo) 
        r[1:3] = pos[1:3]
        u = reshape(pos[4:end],3,3)
        
        w = gramschmidt(u')
        push!(norm1, norm(w[:,1]))
        push!(norm2,norm(w[:,2]))
        push!(norm3, norm(w[:,3]))
        w[:,1] = w[:,1]/norm(w[:,1])
        w[:,2] = w[:,2]/norm(w[:,2])
        w[:,3] = w[:,3]/norm(w[:,3])
        r[4:end] = copy(w'[:])

    end

    exp1 = sum(log(norm1))/(nsteps*dtsampling)
    exp2 = sum(log(norm2))/(nsteps*dtsampling)
    exp3 = sum(log(norm3))/(nsteps*dtsampling)

    println("Exponentes de Lyapunov: $exp1, $exp2, $exp3")

    return r, norm1, norm2, norm3, exp1, exp2, exp3

end

function gramschmidt(u::Matrix{Float64})
     w = eye(3)
     w[:,1] = u[:,1];
     v1 = w[:,1]/norm(w[:,1])
     w[:,2] = u[:,2] - dot(u[:,2],v1)*v1;
     v2 = w[:,2]/norm(w[:,2]);
     w[:,3] = (u[:,3] - dot(u[:,3],v2)*v2 - dot(u[:,3],v1)*v1)
    
    return w
end


    
    
    

    
