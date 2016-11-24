##Recall that rho(p,q,z) = exp(-beta*H)*f(z)

"""
      Random initial vector with components (q,p,zeta)
                """
function initialcondition(beta::Float64, Q::Float64)
    sigma = 1./sqrt(beta)  #Standard deviation
    q,p = rand(Normal(0.0, sigma), 2)
    z = rand(Logistic(0.0,Q))
    return [q,p,z]
end




