using HDF5
using PyPlot

"""
    Open files with the same parameters Q and T
    file1 = h5open("../data/potential/filename.hdf5","r")
    file2 = h5open("../data/potential/filename.hdf5","r")
    ...
    filen = h5open("../data/potential/filename.hdf5","r")
"""

potential = "MH"
cd("../data/$potential")
figname = "lyap10000"*potential
expsinitcond = zeros(1,6)
for i in filter(x -> endswith(x,".hdf5"), readdir())
    filename = "$i"
    file = h5open("$filename","r")
    l = read(file, "lyapexpsinitcond")
    expsinitcond = vcat(expsinitcond,l)
end

l = expsinitcond[2:end,:]


try
     mkdir("../../plots")
end

fig = plt[:figure](figsize=(6,8))
fig[:subplots_adjust](hspace=.5)


ax = fig[:add_subplot](311)
ax[:set_xlabel](L"$\lambda_1$",fontsize="18")
ax[:set_ylabel](L"$q$",fontsize="18")
ax[:hist](l[:,1])

ax = fig[:add_subplot](312)
ax[:set_xlabel](L"$\lambda_2$",fontsize="18")
ax[:set_ylabel](L"$p$", fontsize="18")
ax[:hist](l[:,2])

ax = fig[:add_subplot](313)
ax[:set_xlabel](L"$\lambda_3$",fontsize="18")
ax[:set_ylabel](L"$S$",fontsize="18")
ax[:hist](l[:,3])

plt[:savefig]("../../plots/$figname.png")

writedlm("./$figname", l)

println("Figure succesfully generated. See file  ../plots/$figname.png")
