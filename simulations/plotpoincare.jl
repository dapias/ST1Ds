include("../src/poincare.jl")

using PyPlot

try
    mkdir("../plots")
end

println("Type the name of the file (without a format)")
input = string(readline(STDIN))
filename = input[1:end-1]
potentialname = input[end-2:end-1]

file = h5open("../poincaredata/$potentialname/$(filename).hdf5","r")
data = read(file, "tx")
Q = read(attrs(file)["Q"])
T = read(attrs(file)["T"])
potentialname = read(attrs(file)["potential"])

if potentialname == "HO"
    potential(x) = x^2/2.
elseif potentialname == "QP"
    potential(x) = x^4/4.
elseif potentialname == "MH"
    potential(x) = -1/2.*x^2 + 1/4.*x^4
elseif potentialname == "MP"
    potential(x) =  (1.-exp(-0.5*x))^2.
end


beta = 1./T
ps = compute(data, Q, beta, potential)
ps2 = compute(data, Q, beta, potential, -1.)
plot(ps[:,2],ps[:,4], ".")
plot(ps2[:,2],ps2[:,4], ".")
plt[:xlabel](L"q", fontsize = 18)
plt[:ylabel](L"\zeta", fontsize = 18)

plt[:savefig]("../plots/$filename.png")


close(file)

println("Plot $(filename).png succesfully generated. See file in ../plots")
