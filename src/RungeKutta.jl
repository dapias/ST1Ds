"""
    rungeK(xinit, vectorfield, stepsize)
Classical Runge-Kutta method RK4
"""
function rungeK(x_init::Vector{Float64},field::Function,stepsize::Float64)
  state = x_init
  k1 = field(state)
  k2 = field(state + .5*stepsize*k1)
  k3 = field(state + .5*stepsize*k2)
  k4 = field(state + stepsize*k3)
  phi = (1/6.)*k1 + (1/3.)*k2 + (1/3.)*k3 + (1/6.)*k4
  state = state+stepsize*phi
  return state
end
