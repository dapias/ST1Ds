using HDF5



function makelyaparrays(exps::Int64, potential::AbstractString)

    ##Number of simulation files
    cd("$(homedir())/HooverPrize/data/$potential")
    dir = filter(x -> endswith(x,".hdf5"), readdir())
    sims = length(dir)

    ##Initialize arrays lyap1, ...
    for i in 1:exps
        m = Symbol("lyap$i")
        @eval $m = @eval $(zeros(sims))
    end

    ##Extract data of each file
    j = 1
    for i in dir
        filename = "$i"
        file = h5open("$filename", "r");
        sim = read(file);
        for k in 1:exps
            p = Symbol("lyap$k")
            lyap = @eval $p
            lyap[j]  = @eval ($(sim["exp$k"]))
        end
        j += 1
    end
    return lyap1, lyap2, lyap3, sims
end

println("Type the type of the potential (Harmonic oscillator (HO), Mexican Hat (MH), Quartic (Quartic)) ")
input = string(readline(STDIN))
potential = input[1:end-1]

potentiallist = ["HO", "MH", "Quartic"]

while !(potential in potentiallist)
  println("The potential you typed is not in our database. Try one of the following: \n HO, MH or Quartic or check the spelling")
  input = string(readline(STDIN))
  potential = input[1:end-1]
end


exponents = 3
lyap1, lyap2, lyap3, sims  = makelyaparrays(exponents, potential)
lyapexps = [lyap1, lyap2, lyap3]

cd("$(homedir())/HooverPrize/analysis/")

writedlm("lyapexps$potential.txt", lyapexps)

##Open it via readdlm("lyapexps.txt")
