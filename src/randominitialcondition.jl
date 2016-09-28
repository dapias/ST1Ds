##Recall that rho(p,q,S) = exp(-beta*H)*f(S)

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
function boxmuller(beta::Float64) ##For the harmonic oscillator
    sigma = 1./sqrt(beta)
    u1 = rand()
    u2 = rand()

    p = sqrt(-2*log(u1))*cos(2*pi*u2)
    q = sqrt(-2*log(u1))*sin(2*pi*u2)

    p = p*sigma
    q = q*sigma

    q,p
end
    
"""
Random initial vector with components (q,p,zeta)
"""
function initcond(beta::Float64, Q::Float64)
    q, p = boxmuller(beta)
    s = logisticsampling(Q)
    z = Vector{Float64}([q,p,s])
end



    
