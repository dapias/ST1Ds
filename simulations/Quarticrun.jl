using HDF5
include("./Quarticsimulation.jl")

nsimulations = 10000
T = 10.0
Q = 2.0
nsteps = 100000
deltatsampling = 0.05
deltat = 0.0125


try
    mkdir("../data/Quartic/")
end

for i in 1:nsimulations

    exp1,exp2,exp3 = simulation(T,Q,nsteps,deltatsampling,deltat)
    
    filename = randstring()
    file = h5open("../data/Quartic/$filename.hdf5", "w")

    file["exp1"] = exp1
    file["exp2"] = exp2
    file["exp3"] = exp3

    attrs(file)["nsteps"] = nsteps
    attrs(file)["dtsampling"] = deltatsampling
    attrs(file)["dtintegration"] = deltat


    close(file)

    println("Simulation$i: $filename")
end
