include("../src/solver.jl")
include("./simulation.jl")

using HDF5
using YAML

println("Type the type of the potential (Harmonic oscillator (HO), Mexican Hat (MH), Quartic (QP)) ")
input = string(readline(STDIN))
potentialname = input[1:end-1]

potentiallist = ["HO", "MH", "QP"]

while !(potentialname in potentiallist)
  println("The potential you typed is not in our database. Try one of the following: \n HO, MH or QP or check the spelling")
  input = string(readline(STDIN))
  potentialname = input[1:end-1]
end

pot = potential(potentialname)

try
    mkdir("../poincaredata/")
end

try
    mkdir("../poincaredata/$potentialname/")
end


parameters = YAML.load(open("parametersintegration.yaml"))
T = parameters["T"]
Q = parameters["Q"]
nsteps = parameters["nsteps"]
dt = parameters["dt"]
nsimulations = parameters["nsimulations"]

beta = 1./T;
r0 = initcond(beta, Q)
tfinal= nsteps*dt
filename = randstring(4)
filename *= potentialname
filenamei = filename*"1"

file = h5open("../poincaredata/$potentialname/$(filenamei).hdf5", "w")
(t, xsol) = flowode45(DDfield, r0,dt, tfinal, pot, beta, Q)
x = map(v -> v[1], xsol)
y = map(v -> v[2], xsol)
z = map(v -> v[3], xsol);
tx = [t x y z]
file["tx"] = tx
attrs(file)["Q"] = Q
attrs(file)["T"] = T
attrs(file)["potential"] = "$potentialname"
attrs(file)["dt"] = dt
attrs(file)["nsteps"] = nsteps          
r0 = tx[end,:][2:end]
close(file)
println("Simulation 1 done. File $(filenamei).hdf5 ")

for i in 2:nsimulations
    filenamei = filename*"$i"
    file = h5open("../poincaredata/$potentialname/$(filenamei).hdf5", "w")
    (t, xsol) = flowode45(DDfield, r0,dt, tfinal, pot, beta, Q)
    x = map(v -> v[1], xsol)
    y = map(v -> v[2], xsol)
    z = map(v -> v[3], xsol)
    tx = [t x y z]
    file["tx"] = tx
    r0 = tx[end,:][2:end]
    close(file)
    println("Simulation $i done. File $(filenamei).hdf5 ")
end
 
   
 println("Sequence $(filename)-i.hdf5 succesfully generated. See files in ../poincaredata/$(potentialname)")

#To call the file
# file = h5open("data.hdf5", "r")
# data = read(file, "tx")
