"""
    lyapunov_flow(vectorfield, r0, parameters)
Integrates the equation of motion of a given vector field with the integrator passed through the parameters.
"""
function lyapunov_flow(field::Function, r0::Vector{Float64}, p::Parameters)
    method = p.integrator.f
    t = 0.0:p.dt:p.dtsampling
    
    function extendedfield(time::Float64, r::Vector{Float64})
        field(r, p.potential, 1/p.T, p.thermo)
    end
    
    try
        t,pos = method(extendedfield, r0, t, points=:specified)
        return t, pos
    catch
        t,pos = method(extendedfield, r0, t)
        return t,pos
    end

end


"""
Return the lyapunovspectrum of the passed vector field computed under the conditions passed in the parameters
"""
function lyapunov_spectrum(field::Function, r::Vector{Float64}, p::Parameters, burninsteps = 100)
    n = 3 ##Dimension of the dynamical system considered
    w = eye(n)
    norms = zeros(p.nsteps,n)

    for i in 1:burninsteps
        r = lyapunov_flow(field, r, p)[2][end]
    end

    for i in 1:p.nsteps
        r = lyapunov_flow(field, r, p)[2][end]
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
        exps[k] = sum(log(norms[:,k]))/(p.nsteps*p.dtsampling)
    end

    return exps
end

"""
Return the results of the simulation. Basically, the results are of two types: either the lyapunov spectra for a certain number of random initial conditions or the points belonging to the trajectory defined by a random initial condition.
"""
function lyapunov_simulation(p::Parameters)
    beta = 1./p.T
    r = zeros(12)
    init = initialcondition(beta,p.Q)
    r[1:3] = copy(init)
    r[4] = r[8] = r[12] = 1.0  #Entries of the identity matrix in the initial condition
    exps = lyapunov_spectrum(variationalDDfield,r,p)
    return vcat(exps, init)
end

