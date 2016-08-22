using HDF5

exps = 3
sims = 6

##Initialize arrays C_vGaussian, C_vLogistic ...
for i in 1:exps
    m = symbol("lyap$i")
    @eval ($m = zeros(sims))
end




##Extract data of each file
cd("$(homedir())/HooverPrize/data/HO")
j = 1
for i in filter(x -> endswith(x,".hdf5"), readdir()) #Filter the files that start with string temp
    filename = "$i"
    file = h5open("$filename", "r");
    sim = read(file);
    for k in 1:exps
        p = Symbol("lyap$k")
        @eval ($p[j]) = @eval ($(sim["exp$k"]))
    end
    j += 1
end



cd("$(homedir())/HooverPrize/analysis/")
writedlm("lyap1.txt", lyap1)

file = open("lyap1.txt","a")

writedlm(file, lyap2)
writedlm(file, lyap3)

close(file)
