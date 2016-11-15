const zerotol = 1.e-12


"""
    section(trajectory, thermo, beta, potential, plane, direction)
Compute the intersection of a given trajectory (time-position) with the plane `plane=0`, i.e. plane = "p","q" or "z".
The sense of the intersection is given by the  direction (+1. or -1.).
"""

function section(tx::Matrix{Float64},  thermo::Thermostat, beta::Float64, potential::Potential, plane::String, direction = 1.)

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




"""
    This function computes the Poincare sections for the integrated trajectory saved in trajectorydata. It generates a .hdf5 file.

    Example:
    ```
    julia> sections("myfile", Potential("quartic", x -> x^4/4.), Thermostat("logistic", Q, z-> exp(z/Q)/(Q*( 1 +exp(z/Q))^2.)))
    ```
    """

function sections(filename::String, potential::Potential, thermo::Thermostat)
    filename = filename[end-4:end]
    file = h5open("../data/$filename.hdf5")
    tx = read(file["tx"])
    T = read(attrs(file)["T"])
    nsimulations = read(attrs(file)["nsimulations"])
    beta = 1./T
    close(file)

    zs1 = section(tx, thermo, beta, potential, "z")
    zs2 = section(tx, thermo, beta, potential, "z",-1.)
    ps1 = section(tx, thermo, beta, potential, "p")
    ps2 = section(tx, thermo, beta, potential, "p",-1.)

    zs = vcat(zs1, zs2)
    ps = vcat(ps1, ps2)

      
    newfile =  h5open("../data/$filename.hdf5", "w")

 
    newfile["zsection"] = zs
    newfile["psection"] = ps

    attrs(newfile)["T"] = T

    close(newfile)

    println("File $filename.hdf5 succesfully generated. See file in ../data/")

end
