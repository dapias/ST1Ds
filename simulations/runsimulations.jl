println("Type the type of the potential (Harmonic oscillator (HO), Mexican Hat (MH), Quartic (QP)) ")
input = string(readline(STDIN))
potential = input[1:end-1]

potentiallist = ["HO", "MH", "QP"]

while !(potential in potentiallist)
  println("The potential you typed is not in our database. Try one of the following: \n HO, MH or QP or check the spelling")
  input = string(readline(STDIN))
  potential = input[1:end-1]
end

using HDF5
using YAML

include("./$(potential)simulation.jl")

parameters = YAML.load(open("parameterssimulation.yaml"))

T = parameters["T"]
Q = parameters["Q"]
nsimulations = parameters["nsimulations"]
nsteps = parameters["nsteps"]
deltatsampling = parameters["deltatsampling"]
deltat = parameters["deltat"]


try
    mkdir("../data/")
end

try
    mkdir("../data/$potential/")
end

filename = randstring(4)
file = h5open("../data/$potential/$(filename)$(potential).hdf5", "w")


attrs(file)["nsteps"] = nsteps
attrs(file)["dtsampling"] = deltatsampling
attrs(file)["dtintegration"] = deltat
attrs(file)["Q"] = Q
attrs(file)["T"] = T

close(file)

for i in 1:nsimulations

    file = h5open("../data/$potential/$(filename)$(potential).hdf5", "r+")

    init, exp1,exp2,exp3 = simulation(T,Q,nsteps,deltatsampling,deltat)
    
    file["simulation-$i/initialcond"] = init
    file["simulation-$i/exp1"] = exp1
    file["simulation-$i/exp2"] = exp2
    file["simulation-$i/exp3"] = exp3

    println("Simulation$i done")

    close(file)
end

println("File $(filename)$(potential).hdf5 succesfully generated. See file in ../data/$(potential)")
