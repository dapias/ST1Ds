function thermostatstep(r::Vector{Float64}, par::Parameters)
    q,p,z,v = r
    beta = 1./par.T
    psquare = p^2.
    z += par.dt/4.*(psquare - 1/beta)
    v -= par.dt/2.*(friction(z, par.thermo)/beta)
    p = p*exp(par.dt/2.*(friction(z,par.thermo)/beta))
    z += par.dt/4.*(exp(par.dt*friction(z,par.thermo)/beta)*psquare - 1./beta)

    return [q,p,z,v]
end

function verletstep(r::Vector{Float64}, par::Parameters)
    q,p,z,v = r
    p += 0.5*par.dt*force(q,par.potential)
    q += par.dt*p
    p += 0.5*par.dt*force(q,par.potential)

    return [q,p,z,v]
end


function geometricintegration(r0::Vector{Float64}, p::Parameters)
    t = 0.0:p.dt:p.nsteps*p.dt

    ra = thermostatstep(r0,p)
    rb = verletstep(ra,p)
    r = thermostatstep(rb,p)

    pos = [r0, r]

    for i in 2:p.nsteps
        ra = thermostatstep(r,p)
        rb = verletstep(ra,p)
        r = thermostatstep(rb,p)
        push!(pos,r)
    end

    sample = round(Int,p.dtsampling/p.dt)
    t,xsol = t[1:sample:end], pos[1:sample:end]

    x = map(v -> v[1], xsol)
    y = map(v -> v[2], xsol)
    z = map(v -> v[3], xsol)
    v = map(v -> v[4], xsol)
    tx = [t x y z v]
    return tx
    
end



        
   

