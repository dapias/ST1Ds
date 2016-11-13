using ODE

results = "trajectory" #Two choices: lyapunov or trajectory
T = 1.0
Q = 0.1
dtsampling = 0.0125 
dt=0.0125
nsteps = 10000 #If lyapunov is chosen. Total time = nsteps*dt
nsimulations = 2 #If trajectory is chosen. Total time = nsteps*dt*nsimulations. The data is saved in nsimulations files. 
thermo = Thermostat("logistic", Q, z-> exp(z/Q)/(Q*( 1 +exp(z/Q))^2.))
potential = Potential("quartic", x->x^4./4.)
integrator = Integrator("RK45", ode45)

