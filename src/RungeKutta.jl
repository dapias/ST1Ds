function rungeK(x_init,campo,h)
  state = x_init
  k1 = campo(state)
  k2 = campo(state + .5*h*k1)
  k3 = campo(state + .5*h*k2)
  k4 = campo(state + h*k3)
  phi = (1/6.)*k1 + (1/3.)*k2 + (1/3.)*k3 + (1/6.)*k4
  state = state+h*phi
  return state
end
