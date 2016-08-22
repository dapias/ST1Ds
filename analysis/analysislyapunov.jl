using HDF5

function makelyaparrays(exps::Int64)

    ##Number of simulation files
    cd("$(homedir())/HooverPrize/data/HO")
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

exponents = 3
lyap1, lyap2, lyap3, sims  = makelyaparrays(exponents)
lyapexps = [lyap1, lyap2, lyap3]

cd("$(homedir())/HooverPrize/analysis/")

writedlm("lyapexps.txt", lyapexps)

##Open it via readdlm("lyapexps.txt")
