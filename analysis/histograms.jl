using HDF5

println("Type the filename")
input = string(readline(STDIN))
filename = input[1:end-1]
potentialname=filename[end-1:end]

file = h5open("../poincaredata/sectionsandtrajectories/$filename.hdf5","r")

file2 = h5open("../poincaredata/$potentialname/$(filename)1.hdf5","r")

data = read(file["trajectory"])

q = data[:,2][1:10:end]
p = data[:,3][1:10:end]
z = data[:,4][1:10:end]
Q = read(attrs(file2)["Q"])
T = read(attrs(file2)["T"])
beta = 1./T

    if potentialname == "HO"
        potential(x) = x^2/2.
	rhoq =  sqrt(beta/(2.*pi))*exp(-beta*Float64[potential(i) for i in q])	    
    elseif potentialname == "QP"
        potential(x) = x^4/4.
	rhoq =  1/2.563693352*exp(-beta*Float64[potential(i) for i in q])	    
    elseif potentialname == "MH"
        potential(x) = -1/2.*x^2 + 1/4.*x^4
	rhoq =  1./3.905137170*exp(-beta*Float64[potential(i) for i in q])	    
    end

    
    rhop = sqrt(beta/(2.*pi))*exp(-beta*p.^2/2.)
    rhoz = exp(z/Q)./(Q*(1.+ exp(z/Q)).^2);

q1 = hcat(q, rhoq)
p1 = hcat(p, rhop)
z1 = hcat(z, rhoz)

results = hcat(q1, p1)
results = hcat(results,z1)

writedlm("../data/hist$filename", results)
