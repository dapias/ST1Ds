include("../src/lyapunovspectra.jl")

function simulation(T::Float64,Q::Float64, nsteps::Int64,deltatsampling::Float64, deltat::Float64)
    
    potential(x::Real) = 1/4.*x^4

    beta = 1./T;

    #Burn-in period
    r = zeros(12)
    init = initcond(beta,Q)
    r[1:3] = copy(init)
    r[4] = r[8] = r[12] = 1.0  #Entries of the identity matrix in the initial condition
    burn_intime = 20.0
 
    r = flowRK(variationalDDfield, r, deltat,burn_intime, potential, beta, Q)  #New initial condition after a transient time
    r[4:end] = [1.,0.,0.,0.,1.,0.,0.,0.,1.] #Entries of the identity matrix in the new initial condition;

    ##Simulation to calculate Lyapunov Exponents
    norm1, norm2, norm3, exp1, exp2, exp3 = lyapunovspectra(variationalDDfield,r,deltat, deltatsampling,nsteps, potential, beta, Q)

    return init, exp1,exp2,exp3

end

