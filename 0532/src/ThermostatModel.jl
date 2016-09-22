module ThermostatModel

export Thermostat, Gaussian, Logistic, Quartic, logrhoextended, friction, difffriction

abstract Thermostat

############Nosé-Hooover ###############

type Gaussian{T} <: Thermostat
    Q::T #Parameter that characterizes the distribution (standard deviation)
    beta::T
end

@doc """Logarithm of the distibution in the extended space (f(w))"""->
function logrhoextended(T::Gaussian, x::Float64)
    -T.beta*x^2.0/(2.*T.Q)
end

@doc """Friction coefficient f'(w)/f(w)"""->
function friction(T::Gaussian, x::Float64)
    return -T.beta*x/T.Q
end

@doc """Derivative of the friction coefficient"""->
function difffriction(T::Gaussian, x::Float64)
    return -T.beta/T.Q
end
####################################################


############Logistic distribution##############

type Logistic{T} <: Thermostat
    Q::T ##Parameter that characterizes the distribution (mean)
    beta::T
end

function logrhoextended(T::Logistic, x::Float64)
    z = x-T.Q
    distribution = exp(z)/((1+exp(z))^2)
    log(distribution)
end

@doc """Friction coefficient f'(w)/f(w)"""->
function friction(T::Logistic, x::Float64)
    z = x-T.Q
    (1 - exp(z))/(1+exp(z))
end

@doc """Derivative of the friction coefficient"""->
function difffriction(T::Logistic, x::Float64)
    z = x-T.Q
    return (-2.0*exp(z))/(1.0+exp(z))^2.
end


##############################################

############Quartic Distribution (Fukuda) ###############
type Quartic{T} <: Thermostat
    Q::T ##Parameter that characterizes the distribution
    beta::T
end


# # @doc """Logarithm of the distibution in the extended space (f(w))"""->
function logrhoextended(T::Quartic, x::Float64)
    -T.Q*x^4.
end

@doc """Friction coefficient f'(w)/f(w)"""->
function friction(T::Quartic, x::Float64)
    return -4.0*T.Q*x^3.0
end

@doc """Derivative of the friction coefficient"""->
function difffriction(T::Quartic, x::Float64)
    return -12.0*T.Q*x^2.0
end

end
