include("./DDfield.jl")

using ODE
using HDF5

const zerotol = 1.e-12


"""
    computepsection(trajectory, Q, beta, potential, direction)

Compute the intersection of a given trajectory (time-position) with the plane `p=0`. The sense of the intersection is given by the integer direction (+1 or -1).
"""
function computepsection(tx::Matrix{Float64}, Q::Float64, beta::Float64, potential::Function, direction = 1.)

    function DDextended(time, r::Vector{Float64})
        DDfield(r, potential, beta, Q)
    end

    xsol = tx[:,2:4]  #Spatial part of the trajectory
    e_z = [0.,0.,1.]  #Unit zeta vector
    e_q = [1.,0.,0.]  #Unit q vector

    nhat = cross(e_z, e_q)   #Unit p vector


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
    computezsection(trajectory, Q, beta, potential, direction)

Compute the intersection of a given trajectory (time-position) with the plane `z=0`. The sense of the intersection is given by the integer direction (+1 or -1).
"""

function computezetasection(tx, Q, beta, potential, direction = 1.)

    function DDextended(time, r::Vector{Float64})
        DDfield(r, potential, beta, Q)
    end

    xsol = tx[:,2:4]
    e_p = [0.,1.,0.]
    e_q = [1.,0.,0.]

    nhat = cross(e_p, e_q)

    function projection(x)
            #Projection of a point x over the Poincare section
            dot(x, nhat)*direction
    end

    projarray = [projection(xsol[i,:]) for i in 1:length(xsol[:,1])]

    ps = zeros(4)'

    for i in 1:(length(projarray)-1)
        if  projarray[i] < 0. && projarray[i+1] > 0.
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

                
                
