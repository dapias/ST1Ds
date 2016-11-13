using HDF5

include("../src/DDfield.jl")
include("newparameters.jl")

try
    mkdir("../lyapunovdata/")
end

try
    mkdir("../lyapunovdata/$(potential.name)/")
end


try
    mkdir("../trajectorydata/")
end

try
    mkdir("../trajectorydata/$(potential.name)/")
end


parameters = Parameters(results, T, Q, dtsampling, dt, nsimulations, nsteps, thermo, potential, integrator)

function writeattributes(file, p::Parameters)
    attrs(file)["Integrator"] = p.integrator.name
    attrs(file)["Thermostat"] = p.thermo.name
    attrs(file)["nsteps"] = p.nsteps
    attrs(file)["dtsampling"] = p.dtsampling
    attrs(file)["dtintegration"] = p.dt
    attrs(file)["Q"] = p.Q
    attrs(file)["T"] = p.T
    attrs(file)["nsimulations"] = p.nsimulations
end
    
function run(p::Parameters)
    if p.results == "lyapunov"
        
        filename = randstring(5)    
        file = h5open("../$(p.results)data/$(potential.name)/$(filename).hdf5", "w")
        writeattributes(file,p)
        close(file)
        
        for i in 1:nsimulations
            file = h5open("../$(p.results)data/$(potential.name)/$(filename).hdf5", "r+")
            
            init, exp1,exp2,exp3 = simulation(p)
            file["simulation-$i/initialcond"] = init
            file["simulation-$i/exp1"] = exp1
            file["simulation-$i/exp2"] = exp2
            file["simulation-$i/exp3"] = exp3
            println("Simulation$i done")
            close(file)
        end
        
        println("File $(filename) succesfully generated. See file in ../$(p.results)data/$(potential.name)")

    elseif p.results == "trajectory"

        filename = randstring(5)
        filenamei = filename*"1"
        file = h5open("../$(p.results)data/$(potential.name)/$(filenamei).hdf5", "w")
        writeattributes(file,p)
        tx = simulation(p)
        file["tx"] = tx
        r0 = tx[end,:][2:end]
        close(file)
        println("Simulation 1 done. File $(filenamei).hdf5 ")

        for i in 2:nsimulations
            filenamei = filename*"$i"
            file = h5open("../$(p.results)data/$(potential.name)/$(filenamei).hdf5", "w")
            (t, xsol) = flow(DDfield, r0,p)
            x = map(v -> v[1], xsol)
            y = map(v -> v[2], xsol)
            z = map(v -> v[3], xsol)
            tx = [t x y z]
            file["tx"] = tx
            r0 = tx[end,:][2:end]
            close(file)
            println("Simulation $i done. File $(filenamei).hdf5 ")
        end
        println("Sequence $(filename)-i.hdf5 succesfully generated. See files in ../$(p.results)data/$(potential.name)")
    end
    
end

run(parameters)
