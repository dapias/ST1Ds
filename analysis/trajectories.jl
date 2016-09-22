include("../src/poincare.jl")
using HDF5

println("Type the fileseriesname:")
input = string(readline(STDIN))
fileseriesname = input[1:end-1]

#fileseriesname = "wnvUMH"
potentialfile = fileseriesname[end-1:end]
nsimulations = 1000

i = 1
filename = fileseriesname*"$i"
file = h5open("../poincaredata/$potentialfile/$filename.hdf5")
tx = read(file["tx"])
Q = read(attrs(file)["Q"])
T = read(attrs(file)["T"])
beta = 1./T
potentialname = read(attrs(file)["potential"])
close(file)

if potentialname == "HO"
    potential(x) = x^2/2.
elseif potentialname == "QP"
    potential(x) = x^4/4.
elseif potentialname == "MH"
    potential(x) = -1/2.*x^2 + 1/4.*x^4
end

println("comence bien")
zs1 = computezetasection(tx, Q, beta, potential)
zs2 = computezetasection(tx, Q, beta, potential, -1.)
ps1 = compute(tx, Q, beta, potential)
ps2 = compute(tx, Q, beta, potential, -1.)
println("Section 1 done")

for i in 2:nsimulations
    filename = fileseriesname*"$i"
    file = h5open("../poincaredata/$potentialfile/$filename.hdf5",  "r")
    data = read(file["tx"])
    data = data[2:end,:]
    s1 = computezetasection(data, Q, beta, potential)
    s2 = computezetasection(data, Q, beta, potential, -1.)
    tx = vcat(tx,data)
    zs1 = vcat(zs1, s1)
    zs2 = vcat(zs2, s2)
    s1 = compute(data, Q, beta, potential)
    s2 = compute(data, Q, beta, potential, -1.)
    ps1 = vcat(ps1, s1)
    ps2 = vcat(ps2, s2)
    close(file)
    println("Section $i done")
end

dt = tx[:,1][2]
t = [i*dt for i in 0:length(tx[:,1])-1];
tx[:,1] = t

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





