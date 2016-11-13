include("../src/computesections.jl")
using HDF5

"""
This file analyzes the series of files obtained by integrated a certain initial condition that was saved in trajectorydata. It generates a .hdf5 file with the Poincare sections and the whole trajectory of the integrated solution. The data is grouped in different datasets of the same file.
"""

##Modify this part accordingly
fileseriesname = "oFTyc"
potentialname = "quartic"
potential = Potential(potentialname, x-> x^4/4.)
i = 1
filename = fileseriesname*"$i"
file = h5open("../trajectorydata/$potentialname/$filename.hdf5")
##########

tx = read(file["tx"])
Q = read(attrs(file)["Q"])
T = read(attrs(file)["T"])
nsimulations = read(attrs(file)["nsimulations"])
beta = 1./T
close(file)

thermo = Thermostat("logistic", Q, z-> exp(z/Q)/(Q*( 1 +exp(z/Q))^2.)) 
zs1 = computesection(tx, thermo, beta, potential, "z")
zs2 = computesection(tx, thermo, beta, potential, "z",-1.)
ps1 = computesection(tx, thermo, beta, potential, "p")
ps2 = computesection(tx, thermo, beta, potential, "p",-1.)
println("Section 1 done")

for i in 2:nsimulations
    filename = fileseriesname*"$i"
    file = h5open("../trajectorydata/$potentialname/$filename.hdf5",  "r")
    data = read(file["tx"])
    data = data[2:end,:]
    s1 = computesection(data, thermo, beta, potential, "z")
    s2 = computesection(data, thermo, beta, potential, "z", -1.)
    tx = vcat(tx,data)
    zs1 = vcat(zs1, s1)
    zs2 = vcat(zs2, s2)
    s1 = computesection(data, thermo, beta, potential, "p")
    s2 = computesection(data, thermo, beta, potential, "p",-1.)
    ps1 = vcat(ps1, s1)
    ps2 = vcat(ps2, s2)
    close(file)
    println("Section $i done")
end

dt = tx[:,1][2]
t = [i*dt for i in 0:length(tx[:,1])-1];
tx[:,1] = t

try
    mkdir("../poincaredata")
end

try
    mkdir("../poincaredata/sectionsandtrajectories/")
end

newfile =  h5open("../poincaredata/sectionsandtrajectories/$fileseriesname.hdf5", "w")

newfile["zsection1"] = zs1
newfile["zsection2"] = zs2
newfile["psection1"] = ps1
newfile["psection2"] = ps2
newfile["trajectory"] = tx  

close(newfile)

println("File $fileseriesname.hdf5 succesfully generated. See file in ../poincaredata/sectionsandtrajectories")





