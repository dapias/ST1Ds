include("poincare.jl")

using PyPlot

try
    mkdir("../plots")
end

println("Type the name of the file (without a format)")
input = string(readline(STDIN))
filename = input[1:end-1]

file = h5open("$(filename).hdf5","r")
data = read(file, "tx")
Q = read(attrs(file)["Q"])
beta = read(attrs(file)["beta"])
potentialname = read(attrs(file)["potential"])

if potentialname == "HO"
    potential(x) = x^2/2.
elseif potentialname == "Quartic"
    potential(x) = x^4/4.
elseif potentialname == "MH"
    potential(x) = -1/2.*x^2 + 1/4.*x^4
end



ps = compute(data, Q, beta, potential)
plot(ps[:,2],ps[:,4], ".")
plt[:xlabel](L"q", fontsize = 18)
plt[:ylabel](L"\zeta", fontsize = 18)

plt[:savefig]("../plots/$filename.png")


close(file)

println("Plot $(filename).png succesfully generated. See file in ../plots")
