using HDF5

println("Type the filename")
input = string(readline(STDIN))
filename = input[1:end-1]
potentialname=filename[end-1:end]

file = h5open("../poincaredata/sectionsandtrajectories/$filename.hdf5","r")
file2 = h5open("../poincaredata/$potentialname/$(filename)1.hdf5","r")

zs1 = read(file["zsection1"])
zs2 = read(file["zsection2"])
ps1 = read(file["psection1"])
ps2 = read(file["psection2"])

Q = read(attrs(file2)["Q"])
T = read(attrs(file2)["T"])
beta = 1./T

section1 = hcat(zs1[:,2], zs1[:,3])
section2 = hcat(zs2[:,2], zs2[:,3])

zsection = vcat(section1,section2)

section1 = hcat(ps1[:,2], ps1[:,4])
section2 = hcat(ps2[:,2],ps2[:,4])

psection = vcat(section1, section2)

writedlm("./zsection$filename", zsection)
writedlm("./psection$filename", psection)

close(file)


