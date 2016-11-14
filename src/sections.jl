include("../src/computesections.jl")


"""
    This function computes the Poincare sections for the integrated trajectory saved in trajectorydata. It generates a .hdf5 file.

    Example:
    ```
    julia> sectionsandtrajectory("myfile", Potential("quartic", x -> x^4/4.), Thermostat("logistic", Q, z-> exp(z/Q)/(Q*( 1 +exp(z/Q))^2.)))
    ```
    """

function sections(filename::String, potential::Potential, thermo::Thermostat)
    file = h5open("../trajectorydata/$(potential.name)/$filename.hdf5")
    tx = read(file["tx"])
    T = read(attrs(file)["T"])
    nsimulations = read(attrs(file)["nsimulations"])
    beta = 1./T
    close(file)

    zs1 = computesection(tx, thermo, beta, potential, "z")
    zs2 = computesection(tx, thermo, beta, potential, "z",-1.)
    ps1 = computesection(tx, thermo, beta, potential, "p")
    ps2 = computesection(tx, thermo, beta, potential, "p",-1.)

    zs = vcat(zs1, zs2)
    ps = vcat(ps1, ps2)

    try
        mkdir("../poincaredata/")
    end

    try
        mkdir("../poincaredata/$(potential.name)")
    end
    
    newfile =  h5open("../poincaredata/$(potential.name)/$filename.hdf5", "w")

 
    newfile["zsection"] = zs
    newfile["psection"] = ps

    attrs(newfile)["T"] = T

    close(newfile)

    println("File $filename.hdf5 succesfully generated. See file in ../poincaredata/$(potential.name)")

end
