include("DDfield.jl")

using ODE
using HDF5

zerotol = 1.e-12



function compute(tx, Q, beta, potential, direction = 1.)

    function DDextended(time, r::Vector{Float64})
        DDfield(r, potential, beta, Q)
    end

#    direction = 1.
    xsol = tx[:,2:4]
    e_z = [0.,0.,1.]
    e_2 = [1.,0.,0.]

    nhat = cross(e_z, e_2)

    cond = [dot(xsol[i,:], nhat)*direction for i in 1:length(xsol[:,1])]
    
    function U(x)
        condition = dot(x, nhat)*direction
    end

    ps = zeros(4)
    ps = ps'

    for i in 1:(length(cond)-1)
        if  cond[i] < 0. && cond[i+1] > 0.
            tpsect = tx[i,:][1]
            xcandidate = xsol[i,:]
            deltat = (tx[i+1,:][1] - tx[i,:][1])/2
            psectcond = U(xcandidate)

            while abs(psectcond) > zerotol
                tint = [0.0; deltat]
                (tint, xint) = ode45(DDextended,xcandidate, tint, points = :specified)
                xcandidate = xint[2,:][]
                psectcond = U(xcandidate)
                tpsect += deltat

                if psectcond >0.
                    tpsect -= deltat
                    xcandidate = xint[1,:][]
                    psectcond = U(xcandidate)
                    deltat = deltat/2.
                end
            end

            pselement = [tpsect; xcandidate]
            ps = [ps; pselement']
        end
    end

    ps = ps[2:end,1:end]
end

function computezetasection(tx, Q, beta, potential, direction = 1.)

    function DDextended(time, r::Vector{Float64})
        DDfield(r, potential, beta, Q)
    end

    xsol = tx[:,2:4]
    e_y = [0.,1.,0.]
    e_2 = [1.,0.,0.]

    nhat = cross(e_y, e_2)

    cond = [dot(xsol[i,:], nhat)*direction for i in 1:length(xsol[:,1])]
    
    function U(x)
        condition = dot(x, nhat)*direction
    end

    ps = zeros(4)
    ps = ps'

    for i in 1:(length(cond)-1)
        if  cond[i] < 0. && cond[i+1] > 0.
            tpsect = tx[i,:][1]
            xcandidate = xsol[i,:]
            deltat = (tx[i+1,:][1] - tx[i,:][1])/2
            psectcond = U(xcandidate)

            while abs(psectcond) > zerotol
                tint = [0.0; deltat]
                (tint, xint) = ode45(DDextended,xcandidate, tint, points = :specified)
                xcandidate = xint[2,:][]
                psectcond = U(xcandidate)
                tpsect += deltat

                if psectcond >0.
                    tpsect -= deltat
                    xcandidate = xint[1,:][]
                    psectcond = U(xcandidate)
                    deltat = deltat/2.
                end
            end

            pselement = [tpsect; xcandidate]
            ps = [ps; pselement']
        end
    end

    ps = ps[2:end,1:end]
end

                
                
