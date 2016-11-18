"""
    This file generates an array saved as a plain text that contains the integrated trajectory together with each marginal theoretical distribution. The normalization factor associated with the distribution in `q` is passed as an argument. 
    Example:
    ```
    julia> results = marginaldistributions(filename, Potential("quartic", x->x^4/4.), Thermostat("logistic", A, x -> exp(x/A)/(A*(1.+ exp(x/A)).^2)), normalizationfactor)
    ```
    """

function marginaldistributions(data::Matrix{Float64},p::Parameters, normalizationfactor::Float64)

    beta = 1./p.T
    

    rhoq = 1./normalizationfactor*exp(-beta*Float64[p.potential.f(i) for i in data[:,2]])
    rhop = sqrt(beta/(2.*pi))*exp(-beta*data[:,3].^2/2.)
    rhoz = [p.thermo.distribution(i) for i in data[:,4]]
    q1 = hcat(data[:,2], rhoq)
    p1 = hcat(data[:,3], rhop)
    z1 = hcat(data[:,4], rhoz)
    results = hcat(q1, p1)
    results = hcat(results,z1)
    
end

