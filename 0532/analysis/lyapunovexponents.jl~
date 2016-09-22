using HDF5



function makelyaparrays(exps::Int64, potential::AbstractString, filename::AbstractString)

    ##Number of simulation files
    #cd("$(homedir())/lyapunov/data/$potential")
    #dir = filter(x -> endswith(x,".hdf5"), readdir())

    file = h5open("../data/$potential/$filename.hdf5","r+")
    sims = length(names(file))

    initial = Matrix{Float64}(sims,3)

    ##Initialize arrays lyap1, ...
    for i in 1:exps
        m = Symbol("lyap$i")
        @eval $m = @eval $(zeros(sims))
    end

    ##Extract data of each file
    j = 1
    for i in 1:sims
       sim = read(file, "simulation-$i")
        for k in 1:exps
            p = Symbol("lyap$k")
            lyap = @eval $p
            lyap[j]  = @eval ($(sim["exp$k"]))
        end
        initial[j,:] = sim["initialcond"]
        j += 1
    end

    lyapexps = [lyap1 lyap2 lyap3]
    results= hcat(lyapexps, initial)
    file["lyapexpsinitcond"] = results

    close(file)

#    return initial, lyap1, lyap2, lyap3

end

println("Type the type of the potential (Harmonic oscillator (HO), Mexican Hat (MH), Quartic (QP)) ")
input = string(readline(STDIN))
potential = input[1:end-1]

potentiallist = ["HO", "MH", "QP"]

while !(potential in potentiallist)
  println("The potential you typed is not in our database. Try one of the following: \n HO, MH or QP or check the spelling")
  input = string(readline(STDIN))
  potential = input[1:end-1]
end

println("Type the name of the .hdf5 file (without the extension)")
input = string(readline(STDIN))
filename = input[1:end-1]

exponents = 3
makelyaparrays(exponents, potential, filename)
