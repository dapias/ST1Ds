include("../src/solver.jl")

using HDF5
using YAML

println("Type the type of the potential (Harmonic oscillator (HO), Mexican Hat (MH), Morse Potential (MP) or Quartic (QP)) ")
input = string(readline(STDIN))
potentialname = input[1:end-1]

potentiallist = ["HO", "MH", "QP", "MP"]

while !(potentialname in potentiallist)
  println("The potential you typed is not in our database. Try one of the following: \n HO, MH, MP or QP or check the spelling")
  input = string(readline(STDIN))
  potentialname = input[1:end-1]
end

if potentialname == "HO"
    potential(x) = x^2/2.
elseif potentialname == "QP"
    potential(x) = x^4/4.
elseif potentialname == "MH"
    potential(x) = -1/2.*x^2 + 1/4.*x^4
elseif potentialname == "MP"
    potential(x) =  (1.-exp(-0.5*(x-1.)))^2.
end
    

parameters = YAML.load(open("parametersintegration.yaml"))
T = parameters["T"]
Q = parameters["Q"]
nsteps = parameters["nsteps"]
dt = parameters["dt"]

beta = 1./T;
r0 = initcond(beta, Q)
tfinal= nsteps*dt

(t, xsol) = flowode45(DDfield, r0,dt, tfinal, potential, beta, Q)

x = map(v -> v[1], xsol)
y = map(v -> v[2], xsol)
z = map(v -> v[3], xsol);

tx = [t x y z]

try
    mkdir("../poincaredata/")
end

try
    mkdir("../poincaredata/$potentialname/")
end

filename = randstring(4)
file = h5open("../poincaredata/$potentialname/$(filename)$(potentialname).hdf5", "w")

file["tx"] = tx
attrs(file)["Q"] = Q
attrs(file)["T"] = T
attrs(file)["potential"] = "$potentialname"
attrs(file)["dt"] = dt
attrs(file)["nsteps"] = nsteps

            
close(file)

println("File $(filename)$(potentialname).hdf5 succesfully generated. See file in ../poincaredata/$(potentialname)")

#To call the file
# file = h5open("data.hdf5", "r")
# data = read(file, "tx")
