##Recall that rho(p,q,S) = exp(-beta*H)*f(S)

using Distributions
import Distributions.Normal
"""
Logistic distribution with mean 0 sampled by the method of inverse transform sampling
"""
function logisticsampling(Q::Float64)
    u = rand()
    s = -log(1./u - 1.)*Q 
end

"""
Box Muller method to sampling a bivariate gaussian distribution (corresponds to a harmonic oscillator
in the canonical ensemble)
"""
function boxmuller(beta::Float64) 
    sigma = 1./sqrt(beta)
    q,p = rand(Normal(0.0, sigma), 2)
end
    
"""
Random initial vector with components (q,p,zeta)
"""
function initcond(beta::Float64, Q::Float64)
    q, p = boxmuller(beta)
    z = logisticsampling(Q)
    s = Vector{Float64}([q,p,z])
end



    
