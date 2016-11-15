"""
    This function generates an array saved as a plain text that contains the integrated trajectory together with each marginal theoretical distribution. The normalization factor associated with the distribution in `q` is passed as an argument. 
    Example:
    ```
    julia> results = marginaldistributions(filename, Potential("quartic", x->x^4/4.), Thermostat("logistic", A, x -> exp(x/A)/(A*(1.+ exp(x/A)).^2)), normalizationfactor)
    ```
    """

function marginaldistributions(filename::String, potential::Potential, thermo::Thermostat, normalizationfactor::Float64)
    filename = filename[end-4:end]
    file = h5open("../data/$filename.hdf5","r")

    data = read(file["tx"])
    q = data[:,2][1:10:end]
    p = data[:,3][1:10:end]
    z = data[:,4][1:10:end]
    T = read(attrs(file)["T"])
    beta = 1./T

    rhoq = 1./normalizationfactor*exp(-beta*Float64[potential.f(i) for i in q])
    rhop = sqrt(beta/(2.*pi))*exp(-beta*p.^2/2.)
    rhoz = [thermo.distribution(i) for i in z]
    q1 = hcat(q, rhoq)
    p1 = hcat(p, rhop)
    z1 = hcat(z, rhoz)
    results = hcat(q1, p1)
    results = hcat(results,z1)


    writedlm("../data/hist$filename", results)

    println("File hist$filename succesfully generated. See file in ../data/")
    
end

