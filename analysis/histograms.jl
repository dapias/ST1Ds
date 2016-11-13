"""
This file generates an array saved as a plain text that contains the integrated trajectory together with each marginal theoretical distribution
"""

using HDF5

###Modify this part accordingly
filename = "oFTyc"
potentialname = "quartic"
#########

file = h5open("../poincaredata/sectionsandtrajectories/$filename.hdf5","r")
file2 = h5open("../trajectorydata/$potentialname/$(filename)1.hdf5","r")

data = read(file["trajectory"])
q = data[:,2][1:10:end]
p = data[:,3][1:10:end]
z = data[:,4][1:10:end]
Q = read(attrs(file2)["Q"])
T = read(attrs(file2)["T"])
beta = 1./T

####Modify this part in accordance to the used potential
if potentialname == "harmonic"
    potential(x) = x^2/2.
    rhoq =  sqrt(beta/(2.*pi))*exp(-beta*Float64[potential(i) for i in q])	    
elseif potentialname == "quartic"
    potential(x) = x^4/4.
    rhoq =  1/2.563693352*exp(-beta*Float64[potential(i) for i in q])	    
elseif potentialname == "mexican"
    potential(x) = -1/2.*x^2 + 1/4.*x^4
    rhoq =  1./3.905137170*exp(-beta*Float64[potential(i) for i in q])	    
end
#########################3

rhop = sqrt(beta/(2.*pi))*exp(-beta*p.^2/2.)
rhoz = exp(z/Q)./(Q*(1.+ exp(z/Q)).^2)
q1 = hcat(q, rhoq)
p1 = hcat(p, rhop)
z1 = hcat(z, rhoz)
results = hcat(q1, p1)
results = hcat(results,z1)

writedlm("../data/hist$filename", results)

