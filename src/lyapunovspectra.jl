include("./DDfield.jl")
include("./RungeKutta.jl")
include("./randominitialcondition.jl")
include("./gramschmidt.jl")


"""
Solver based on the 4th order Runge-Kutta integrator
"""
function flowRK(field::Function, r0::Vector{Float64},dt::Float64, tfinal::Float64, potential::Function, beta::Float64, Q::Float64)

    t = 0.0:dt:tfinal
    pos = copy(r0)

    function extendedfield(r::Vector{Float64})
        field(r, potential, beta, Q)
    end

    N = length(t) - 1
    for i in 1:N
        pos = rungeK(pos, extendedfield, dt)
    end

    return pos
end

"""
Calculates the whole lyapunov spectrum of the given (variational) field
"""
function lyapunovspectra(field::Function, r::Vector{Float64}, dt::Float64, dtsampling::Float64, nsteps::Int64, potential::Function, beta::Float64, Q::Float64)
    w = eye(3)
    norm1 = zeros(nsteps)
    norm2 = zeros(nsteps)
    norm3 = zeros(nsteps)
#    phasespace = zeros(nsteps,3) #If want to store the trajectory of the point

    for i in 1:nsteps
 #       phasespace[i,:] = r[1:3]
        r = flowRK(field, r, dt, dtsampling, potential, beta, Q) 
#        r[1:3] = pos[1:3]
        
        u = reshape(r[4:end],3,3)
        
        w = gramschmidt(u')
        norm1[i] = norm(w[:,1])
        norm2[i] = norm(w[:,2])
        norm3[i] = norm(w[:,3])

        w[:,1] = w[:,1]/norm(w[:,1])
        w[:,2] = w[:,2]/norm(w[:,2])
        w[:,3] = w[:,3]/norm(w[:,3])
        r[4:end] = copy(w'[:])

    end

    exp1 = sum(log(norm1))/(nsteps*dtsampling)
    exp2 = sum(log(norm2))/(nsteps*dtsampling)
    exp3 = sum(log(norm3))/(nsteps*dtsampling)

#    println("Exponentes de Lyapunov: $exp1, $exp2, $exp3")

    return norm1, norm2, norm3, exp1, exp2, exp3

end



    
    
    

    
