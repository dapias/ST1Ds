const zerotol = 1.e-12


"""
        section(trajectory, thermo, beta, potential, plane, direction)
    Compute the intersection of a given trajectory (time-position) with the plane `plane=0`, i.e. plane = "p","q" or "z".
    The sense of the intersection is given by the  direction (+1. or -1.).
    """

function section(tx::Matrix{Float64},  p::Parameters, plane::String)

    beta = 1/p.T
    
    function DDextended(time, r::Vector{Float64})
        DDfield(r, p.potential, beta, p.thermo)
    end
    xsol = tx[:,2:4]  #Spatial part of the trajectory

    planes = Dict("p"=>[0.,1.,0.], "z"=>[0.,0.,1.],"q"=>[1.,0.,0.])

    nhat = planes[plane]  #Unit vector normal to the chosen plane

    sections = zeros(4)'
    
    for direction in [1.0,-1.0]
        
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
        sections = vcat(sections,ps)
    end

    return sections[2:end,:]
end
