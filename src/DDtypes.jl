module DDtypes

export Thermostat, Potential, Integrator, Parameters

type Thermostat{F<:Function}
    name::String
    Q::Float64
    distribution::F
end

type Potential{F<:Function}
    name::String
    f::F
end

type Integrator{F<:Function}
    name::String
    f::F
end

type Parameters
    results::String
    T::Float64
    Q::Float64
    dtsampling::Float64
    dt::Float64
    nsimulations::Int64
    nsteps::Int64
    thermo::Thermostat
    potential::Potential
    integrator::Integrator
end

end
