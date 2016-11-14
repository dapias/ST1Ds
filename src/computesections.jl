const zerotol = 1.e-12


"""
    computesection(trajectory, thermo, beta, potential, plane, direction)
Compute the intersection of a given trajectory (time-position) with the plane `plane=0`, i.e. plane = "p","q" or "z".
The sense of the intersection is given by the  direction (+1. or -1.).
"""

function computesection(tx::Matrix{Float64},  thermo::Thermostat, beta::Float64, potential::Potential, plane::String, direction = 1.)

  function DDextended(time, r::Vector{Float64})
    DDfield(r, potential, beta, thermo)
  end
  xsol = tx[:,2:4]  #Spatial part of the trajectory

  if plane == "p"
    e_z = [0.,0.,1.]
    e_q = [1.,0.,0.]
    nhat = cross(e_z, e_q)
  elseif plane == "z"
    e_p = [0.,1.,0.]
    e_q = [1.,0.,0.]
    nhat = cross(e_p, e_q)
  elseif plane == "q"
    e_z = [0.,0.,1.]
    e_p = [0.,1.,0.]
    nhat = cross(e_z,e_p)
  end


  function projection(x)
    #Projection of a point x over the Poincare section
    dot(x, nhat)*direction
  end

  projarray = [projection(xsol[i,:]) for i in 1:length(xsol[:,1])]
  ps = zeros(4)'

  for i in 1:(length(projarray)-1)
    if  projarray[i] < 0. && projarray[i+1] > 0. #Condition of intersection
      tpsect = tx[i,:][1]
      xcandidate = xsol[i,:]
      deltat = (tx[i+1,:][1] - tx[i,:][1])/2
      psectcond = projection(xcandidate)

      while abs(psectcond) > zerotol
        tint = [0.0; deltat]
        (tint, xint) = ode45(DDextended,xcandidate, tint, points = :specified)
        xcandidate = xint[2,:][]
        psectcond = projection(xcandidate)
        tpsect += deltat

        if psectcond >0.
          tpsect -= deltat
          xcandidate = xint[1,:][]
          psectcond = projection(xcandidate)
          deltat = deltat/2.
        end
      end

      pselement = [tpsect; xcandidate]
      ps = [ps; pselement']
    end
  end

  ps = ps[2:end,1:end]
end


