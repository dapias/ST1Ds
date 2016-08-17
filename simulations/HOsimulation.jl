include("../src/lyapunovspectra.jl")

potential(x) = x^2/2.

T = 10.0
beta = 1./T;
Q = 2.0; #"mass" ot the thermostat;

#Burn-in period
r = zeros(12)
r[1:3] = initcond(beta, Q) #Initial condition taken from the assumed equilibrium distribution
r[4] = r[8] = r[12] = 1.0  #Entries of the identity matrix in the initial condition
burn_intime = 20.0
deltat = 0.02

r = flowRK(variationalDDfield, r, deltat,burn_intime, potential, beta, Q)  #New initial condition after a transient time
r[4:end] = [1.,0.,0.,0.,1.,0.,0.,0.,1.] #Entries of the identity matrix in the new initial condition;

##Simulation to calculate Lyapunov Exponents
nsteps = 100000
deltatsampling = 0.05
deltat = 0.0125
norm1, norm2, norm3, exp1, exp2, exp3 = lyapunovspectra(variationalDDfield,r,deltat, deltatsampling,nsteps, potential, beta, Q);


