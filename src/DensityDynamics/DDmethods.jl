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
Equations of motion for a general potential coupled to any thermostat
"""
function DDfield(r::Vector{Float64}, potential::Potential, beta::Float64, thermo::Thermostat)

    (q, p, z, v) = r

    f = x-> force(x, potential)

    dq_dt = p
    dp_dt = f(q) + friction(z,thermo)*p/beta
    dz_dt = p^2. - 1.0/beta
    dv_dt = -friction(z,thermo)/beta

    
    [dq_dt, dp_dt, dz_dt, dv_dt]

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




