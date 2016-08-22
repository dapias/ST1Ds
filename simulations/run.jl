println("Type the type of the potential (Harmonic oscillator (HO), Mexican Hat (MH), Quartic (Quartic)) ")
input = string(readline(STDIN))
potential = input[1:end-1]

potentiallist = ["HO", "MH", "Quartic"]

while !(potential in potentiallist)
  println("The potential you typed is not in our database. Try one of the following: \n HO, MH or Quartic or check the spelling")
  input = string(readline(STDIN))
  potential = input[1:end-1]
end

using HDF5

include("./$(potential)simulation.jl")

nsimulations = 5
T = 10.0
Q = 2.0
nsteps = 100000
deltatsampling = 0.05
deltat = 0.0125


try
    mkdir("../data/")
end

try
    mkdir("../data/$potential/")
end

for i in 1:nsimulations

    init, exp1,exp2,exp3 = simulation(T,Q,nsteps,deltatsampling,deltat)
    
    filename = randstring()
    file = h5open("../data/$potential/$filename.hdf5", "w")

    file["initialcond"] = init
    file["exp1"] = exp1
    file["exp2"] = exp2
    file["exp3"] = exp3

    attrs(file)["nsteps"] = nsteps
    attrs(file)["dtsampling"] = deltatsampling
    attrs(file)["dtintegration"] = deltat


    close(file)

    println("Simulation$i: $filename")
end
