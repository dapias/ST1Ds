function geometricmethod(field::Function, r0::Vector{Float64}, p::Parameters)

    
end
    

function trajectoryintegrator(field::Function, r0::Vector{Float64}, p::Parameters)
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
