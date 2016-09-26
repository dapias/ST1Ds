using ForwardDiff
import ForwardDiff.derivative


"""
Equations of motion fot a general potential coupled to the logistic thermostat
"""
function DDfield(r::Vector{Float64}, potential::Function, beta::Float64, Q::Float64)

    (q, p, z) = r

    force(x::Float64) = -derivative(potential,x)
    
    dq_dt = p
    dp_dt = force(q) + (1-exp(z/Q))/(Q*(1+exp(z/Q)))*p/beta
    dz_dt = p^2. - 1.0/beta

    [dq_dt, dp_dt, dz_dt]

end

"""
Returns the derivative of the force given a potential,
i.e. the negative second derivative of the potential
"""
function forcederivative(potential::Function)
    force(x) = -derivative(potential,x)
    fprime(x) = derivative(force,x)

    return fprime

end


"""
    jacobian(r, potential, beta, Q)

Returns the jacobian of the DDfield evaluated at the point `r` for a given `potential` and some values of `beta` and `Q`
"""
function jacobian(r::Vector{Float64}, potential::Function, beta::Float64, Q::Float64)

    (q,p,z) = r
    
    fprime = forcederivative(potential)
    
    J = [0. 1. 0.; fprime(q)  (1-exp(z/Q))/(Q*(1.+exp(z/Q)))/beta p*(-2.0*exp(z/Q))/(Q*(1.0+exp(z/Q)))^2./beta; 0. 2.0*p 0.]

end
    
    
"""
    variationalDDfield(r, potential, beta, Q)

Evaluates the variational field associated to the DDfield at the extended vector `r` (q,p,z  and 9 more entries representing the coordinates of three tangent vectors) for a given `potential` and some values of `beta` and `Q`
"""
function variationalDDfield(r_and_phi::Vector{Float64}, potential::Function, beta::Float64, Q::Float64)

    v = DDfield(r_and_phi[1:3], potential, beta, Q)

    J = jacobian(r_and_phi[1:3], potential, beta, Q)
   
    rmatrix = reshape(r_and_phi[4:end],3,3)

    DPhi = J*rmatrix'
    

    return append!(v, DPhi'[:])
end    
