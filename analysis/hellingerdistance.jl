using HDF5
using StatsBase

println("Type the filename")
input = string(readline(STDIN))
filename = input[1:end-1]
potentialname=filename[end-1:end]

file = h5open("../poincaredata/sectionsandtrajectories/$filename.hdf5","r")
file2 = h5open("../poincaredata/$potentialname/$(filename)1.hdf5","r")

data = read(file["trajectory"])

t = data[:,1]
q = data[:,2]
p = data[:,3]
z = data[:,4];
Q = read(attrs(file2)["Q"])
T = read(attrs(file2)["T"])
beta = 1./T

    if potentialname == "HO"
        potential(x) = x^2/2.
    elseif potentialname == "QP"
        potential(x) = x^4/4.
    elseif potentialname == "MH"
        potential(x) = -1/2.*x^2 + 1/4.*x^4
    end


function hellinger(x,y)
    sqrt(1. - sum(sqrt(x .* y) / sqrt(sum(x) * sum(y))))
end

function jointdistribution(q::Float64,p::Float64, z::Float64, pot::AbstractString)
    if pot == "HO"
        rhoq =  sqrt(beta/(2.*pi))*exp(-beta*potential(q))
    elseif pot == "QP"
        rhoq =  1/2.563693352*exp(-beta*potential(q))
    elseif pot == "MH"
        rhoq =  1./3.905137170*exp(-beta*potential(q))
    end
    
    rhop = sqrt(beta/(2.*pi))*exp(-beta*p.^2/2.)
    rhoz = exp(z/Q)./(Q*(1.+ exp(z/Q)).^2);
    rhoq*rhop*rhoz
end

n = 20
time = zeros(n)
hell = zeros(n)

q = q[1:10:end]
p = p[1:10:end]
z = z[1:10:end]
t = t[1:10:end]
len = length(t)-1 

for i in 1:n
    t1 = t[1:Int(len/n*i)]
    q1 = q[1:Int(len/n*i)]
    p1 = p[1:Int(len/n*i)]
    z1 = z[1:Int(len/n*i)]
    
    time[i] = t1[end]
    
       
    freq_exp = fit(Histogram, (q1,p1,z1), nbins = 100);
    
    edges = freq_exp.edges
    step_x = step(edges[1])
    step_y = step(edges[2])
    step_z = step(edges[3]);
    
    println(i)
    
    hist_exp = freq_exp.weights/(sum(freq_exp.weights*step_x*step_y*step_z));
    
    qrange = edges[1][1]+step_x/2.:step_x:edges[1][end]-step_x/2.
    prange = edges[2][1]+step_y/2.:step_y:edges[2][end]-step_y/2.
    zrange = edges[3][1]+step_z/2.:step_z:edges[3][end]-step_z/2.
    
    hist_theor = Float64[jointdistribution(x,y,z, potentialname) for x in qrange, 
        y in prange, z in zrange];
    
    hell[i] = hellinger(hist_exp, hist_theor)
end

close(file)
close(file2)

hellinger_results = hcat(time, hell)
try
    mkdir("../hellinger/")
end
writedlm("./hell$filename",hellinger_results)




