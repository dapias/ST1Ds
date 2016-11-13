include("./randominitialcondition.jl")
include("./DDtypes.jl")

using ForwardDiff
import ForwardDiff.derivative
using ODE
using DDtypes


function force(x::Float64, potential::Potential)
    f(y) = -derivative(potential.f,y)
    f(x)
end


function friction(z::Float64, thermo::Thermostat)
    g = y -> derivative(thermo.distribution,y)
    g(z)/thermo.distribution(z)
end

function forcederivative(potential::Potential)
    force(x) = -derivative(potential.f,x)
    fprime(x) = derivative(force,x)
end


"""
Equations of motion fot a general potential coupled to any thermostat
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
    jacobian(r, potential, beta, thermo)
Returns the jacobian of the DDfield evaluated at the point `r`
"""
function jacobian(r::Vector{Float64}, potential::Potential, beta::Float64, thermo::Thermostat)

    (q,p,z) = r
    
    fprime = forcederivative(potential)

    distprime(x) = derivative(thermo.distribution,x)
    distprime2(x) = derivative(distprime,x)
    
    
    J = [0. 1. 0.; fprime(q)  friction(z, thermo)/beta p/beta*(distprime2(z)/thermo.distribution(z) - friction(z, thermo)^2.); 0. 2.0*p 0.]

end


"""
Evaluates the variational field associated to the DDfield at the extended vector `r` (q,p,z  and 9 more entries representing the coordinates of three tangent vectors)
"""
function variationalDDfield(r_and_phi::Vector{Float64}, potential::Potential, beta::Float64, thermo::Thermostat)

    v = DDfield(r_and_phi[1:3], potential, beta, thermo)
    J = jacobian(r_and_phi[1:3], potential, beta, thermo)
    rmatrix = reshape(r_and_phi[4:end],3,3)
    DPhi = J*rmatrix'

    return append!(v, DPhi'[:])
end    

"""
    flow(vectorfield, r0, parameters)
Integrates the equation of motion of a given vector field with the integrator passed through the parameters.
"""
function flow(field::Function, r0::Vector{Float64}, p::Parameters)
    method = p.integrator.f  
    if p.results == "lyapunov"
        t = 0.0:p.dt:p.dtsampling
    elseif  p.results == "trajectory"
        t = 0.0:p.dt:p.nsteps*p.dt
    end

    function extendedfield(time::Float64, r::Vector{Float64})
        field(r, p.potential, 1/p.T, p.thermo)
    end

    try
        t,pos = method(extendedfield, r0, t, points=:specified)
    catch
        t,pos = method(extendedfield, r0, t)
    end

    if p.results == "trajectory"
        sample = round(Int,dtsampling/dt)
        t, pos = t[1:sample:end], pos[1:sample:end]
    end

    return t,pos
end


"""
Given a 3x3 matrix  orthogonalize its columns with the Gram-Schmidt procedure
"""
function gramschmidt(u::Matrix{Float64})
    w = eye(3)
    w[:,1] = u[:,1];
    v1 = w[:,1]/norm(w[:,1])
    w[:,2] = u[:,2] - dot(u[:,2],v1)*v1;
    v2 = w[:,2]/norm(w[:,2]);
    w[:,3] = (u[:,3] - dot(u[:,3],v2)*v2 - dot(u[:,3],v1)*v1)
    
    return w
end

"""
Return the lyapunovspectrum of the passed vector field computed under the conditions passed in the parameters
"""
function lyapunovspectrum(field::Function, r::Vector{Float64}, p::Parameters)
    n = 3 ##Dimension of the dynamical system considered
    w = eye(n)
    norms = zeros(p.nsteps,n)

    for i in 1:p.nsteps
        r = flow(field, r, p)[2][end]
        u = reshape(r[4:end],n,n)     
        w = gramschmidt(u')
        for j in 1:n
            norms[i,j] = norm(w[:,j])
        end
        for j in 1:n
            w[:,j] = w[:,j]/norm(w[:,j])
        end
        r[4:end] = copy(w'[:])
    end
    exps = zeros(n)
    for k in 1:n
        exps[k] = sum(log(norms[:,k]))/(nsteps*dtsampling)
    end

    return exps
end

"""
Return the results of the simulation. Basically, the results are of two types: either the lyapunov spectra for a certain number of random initial conditions or the points belonging to the trajectory defined by a random initial condition. 
"""
function simulation(p::Parameters)
    beta = 1./p.T
    if p.results == "lyapunov"
        r = zeros(12)
        init = initcond(beta,p.Q)
        r[1:3] = copy(init)
        r[4] = r[8] = r[12] = 1.0  #Entries of the identity matrix in the initial condition
        burn_intime = 20.0
        r = flow(variationalDDfield, r, p)[end][2]  #New initial condition after a transient time
        r[4:end] = [1.,0.,0.,0.,1.,0.,0.,0.,1.] #Entries of the identity matrix in the new initial condition
        #After burn-in time
        exp1, exp2, exp3 = lyapunovspectrum(variationalDDfield,r,p)

        return init, exp1,exp2,exp3

    elseif p.results == "trajectory"
        r0 = initcond(beta, p.Q)
        (t, xsol) = flow(DDfield, r0, p)
        x = map(v -> v[1], xsol)
        y = map(v -> v[2], xsol)
        z = map(v -> v[3], xsol)
        tx = [t x y z]

        return tx
    end
    
end
