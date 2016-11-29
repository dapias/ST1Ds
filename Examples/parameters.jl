using ODE

T = 1.0
Q = 0.1
dtsampling = 0.0125
dt=0.0025
nsteps = 50000 #If lyapunov_exponent function is chosen, total time = nsteps*dt
nsimulations = 50 #If trajectory function is chosen, total time = nsteps*dt*nsimulations. 
thermo = Thermostat("logistic", z-> exp(z/Q)/(Q*( 1 +exp(z/Q))^2.))
potential = Potential("quartic", x->x^4./4.)
integrator = Integrator("RK4", ode4)
