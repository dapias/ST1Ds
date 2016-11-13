using HDF5
using StatsBase
using Cubature
using PyCall

@pyimport statsmodels.api as sm

###Modify this part accordingly
filename = "oFTyc"
potentialname = "quartic"
#########

file = h5open("../poincaredata/sectionsandtrajectories/$filename.hdf5","r")
file2 = h5open("../trajectorydata/$potentialname/$(filename)1.hdf5","r")

data = read(file["trajectory"])
Q = read(attrs(file2)["Q"])
T = read(attrs(file2)["T"])
beta = 1./T

"""
Returns the theoretical joint distribution given the potential. The logistic distribution in `z` is considered as default
"""
function jointdistribution(x::Vector{Float64}, pot::AbstractString)
    q,p,z = x
    ####Modify this part according to the potential and to the
    #thermostat used
    if pot == "harmonic"
        rhoq =  sqrt(beta/(2.*pi))*exp(-beta*(q^2./2))
    elseif pot == "quartic"
        rhoq =  1/2.563693352*exp(-beta*(q^4./4.))
    elseif pot == "mexican"
        rhoq =  1./3.905137170*exp(-beta*(-q^2./2+q^4/4.))
    end
    rhop = sqrt(beta/(2.*pi))*exp(-beta*p.^2/2.)
    rhoz = exp(z/Q)./(Q*(1.+ exp(z/Q)).^2);
    ############################
    rhoq*rhop*rhoz
end


joint(x::Vector{Float64}) = jointdistribution(x, potentialname)

n = 20  ##Number of points where the hellinger distance is calculated
time = zeros(n)
hell = zeros(n)
error = zeros(n)
data = data[1:100:end,:] ##The 100 factor increases the sampling time by this factor
len = size(data)[1] - 1

for i in 1:n
    data1 = data[1:Int(len/n*i),:]
    time[i] = data1[:,1][end]
    dens = sm.nonparametric[:KDEMultivariate](data = data1[:,2:4], var_type = "ccc")

    function kernel(x::Vector{Float64})
        dens[:pdf](x)[1]
    end

    k_2(x) = kernel(x)^0.5
    j_2(x) = joint(x)^0.5

    hellinger(x) = 2.*(k_2(x) - j_2(x))^2.
    
    q1, q2 = (minimum(data1[:,2]),maximum(data1[:,2]))
    p1, p2 = (minimum(data1[:,3]),maximum(data1[:,3]))
    z1, z2 = (minimum(data1[:,4]),maximum(data1[:,4]))
    hell[i], error[i] = hcubature(hellinger, [q1,p1,z1],[q2,p2,z2], reltol = 1.0e-2)
    println(i)
end

close(file)
close(file2)

results = hcat(time, hell)
hellinger_results = hcat(results, error)
try
    mkdir("../hellinger/")
end
writedlm("../hellinger/integral$filename",hellinger_results)
