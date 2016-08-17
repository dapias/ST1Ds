using HDF5

try
    mkdir("../data/")
    mkdir("../data/HO/")
end

filename = randstring()
file = h5open("../data/HO/$filename.hdf5", "w")

include("./HOsimulation.jl")

file["exp1"] = exp1
file["exp2"] = exp2
file["exp3"] = exp3

attrs(file)["nsteps"] = nsteps
attrs(file)["dtsampling"] = deltatsampling
attrs(file)["dtintegration"] = deltat


close(file)
