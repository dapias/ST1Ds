"""
    trajectory_flow(vectorfield, r0, parameters)
Integrates the equation of motion of a given vector field with the integrator passed through the parameters.
"""
function trajectory_flow(field::Function, r0::Vector{Float64}, p::Parameters)
    method = p.integrator.f
    t = 0.0:p.dt:p.nsteps*p.dt

    function extendedfield(time::Float64, r::Vector{Float64})
        field(r, p.potential, 1/p.T, p.thermo)
    end

    pos = [zeros(3), zeros(3)] #Initialize pos array with the correct type
    
    try
        t,pos = method(extendedfield, r0, t, points=:specified)
    catch
        t,pos = method(extendedfield, r0, t)
    end

    
    sample = round(Int,p.dtsampling/p.dt)
    return t[1:sample:end], pos[1:sample:end]

end



"""
Return the results of the simulation. Basically, the results are of two types: either the lyapunov spectra for a certain number of random initial conditions or the points belonging to the trajectory defined by a random initial condition.
"""
function trajectory_simulation(p::Parameters, r0::Vector{Float64})
    (t, xsol) = trajectory_flow(DDfield, r0, p)
    x = map(v -> v[1], xsol)
    y = map(v -> v[2], xsol)
    z = map(v -> v[3], xsol)
    v = map(v -> v[4], xsol)
    tx = [t x y z v]
    return tx
end

