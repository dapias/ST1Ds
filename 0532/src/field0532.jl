using ForwardDiff
import ForwardDiff.derivative


function field0532(r::Vector{Float64}, potential::Function, beta::Float64, Q::Float64)

    (q, p, z) = r

    force(x::Float64) = -derivative(potential,x)

    
    dq_dt = p
    dp_dt = force(q) -z*(0.05*p + 0.32*p^3*beta)
    dz_dt = 0.05*(p^2.*beta-1.) + 0.32*(p^4*beta^2. - 3.*p^2*beta)

    [dq_dt, dp_dt, dz_dt]

end

function forcederivative(potential::Function)
    force(x) = -derivative(potential,x)
    fprime(x) = derivative(force,x)

    return fprime

end

function jacobian(r_and_phi::Vector{Float64}, potential::Function, beta::Float64, Q::Float64)

    (q, p, z) = r_and_phi[1:3]
    
    fprime = forcederivative(potential)

    
    J = [0. 1. 0.; fprime(q) -z*(0.05 + 3*0.32*p^2.) -(0.05*p + 0.32*p^3); 0.0 (0.1*p + 0.32*(4*p^3. - 6*p)) 0.0]

end
    
    

function variational0532field(r_and_phi::Vector{Float64}, potential::Function, beta::Float64, Q::Float64)

    r = field0532(r_and_phi[1:3], potential, beta, Q)

    J = jacobian(r_and_phi, potential, beta, Q)
   
    rmatrix = reshape(r_and_phi[4:end],3,3)
    DPhi = J*rmatrix'
    

    return append!(r, DPhi'[:])

end    
