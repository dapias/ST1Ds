#include("parameters.jl")

#parameters = Parameters(results, T, Q, dtsampling, dt, nsimulations, nsteps, thermo, potential, integrator)

function writeattributes(file, p::Parameters)
    attrs(file)["Integrator"] = p.integrator.name
    attrs(file)["Thermostat"] = p.thermo.name
    attrs(file)["Potential"] = p.potential.name
    attrs(file)["nsteps"] = p.nsteps
    attrs(file)["dtsampling"] = p.dtsampling
    attrs(file)["dtintegration"] = p.dt
    attrs(file)["Q"] = p.Q
    attrs(file)["T"] = p.T
    attrs(file)["nsimulations"] = p.nsimulations
end

"""
    Returns the results of the simulation in a set of hdf5 files depending on the type of the simulation. If the lyapunov spectra are asked the results are returned in one file where the lyapunov spectrum for each initial condition is saved. If the trajectory is asked a long trajectory will be saved in different files coincident with the number of simulations passed.
    """
function runsimulation(p::Parameters)

    try
        mkdir("../data")
    end
    
    filename = randstring(5)

    function lyapunov()
        file = h5open("../data/lyap$(filename).hdf5", "w")
        writeattributes(file,p)
        close(file)
        
        for i in 1:p.nsimulations
            file = h5open("../data/lyap$(filename).hdf5", "r+")

            init, exp1,exp2,exp3 = simulation(p)
            file["simulation-$i/initialcond"] = init
            file["simulation-$i/exp1"] = exp1
            file["simulation-$i/exp2"] = exp2
            file["simulation-$i/exp3"] = exp3
            println("Simulation$i done")
            close(file)
        end

        println("File lyap$(filename).hdf5 succesfully generated. See file in ../data/")
    end

    function trajectory()
        
        tx = simulation(p)
        r0 = tx[end,:][2:end]
        println("Part 1 done. ")

        for i in 2:p.nsimulations
            (t, xsol) = flow(DDfield, r0,p)
            x = map(v -> v[1], xsol)[2:end]
            y = map(v -> v[2], xsol)[2:end]
            z = map(v -> v[3], xsol)[2:end]
            t += (i-1)*p.nsteps*p.dt
            traj = [t[2:end] x y z]
            tx = vcat(tx, traj)
            r0 = traj[end,:][2:end]
            println("Part $i done.")
        end

        file = h5open("../data/traj$(filename).hdf5", "w")
        file["tx"] = tx
        writeattributes(file,p)
        close(file)

        println("Trajectory traj$(filename).hdf5 succesfully generated. See file in ../data")
    end

    if p.results == "lyapunov"
        lyapunov()
    elseif p.results == "trajectory"
        trajectory()
    end
end

#Executing the main function
#run(parameters)
