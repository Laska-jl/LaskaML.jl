
module hh

const FunctionOrNothing = T where T<:Union{Function, Nothing}
const NumberOrNothing = T where T<:Union{Number, Nothing}
# const ChannelType = NamedTuple{
#         (
#         :id, :g, :Vsteady, :alpha_m, :beta_m, :p, :alpha_h, :beta_h, :q
#         ),
#         Tuple{
#             String,
#             U,
#             T,
#             FunctionOrNothing,
#             FunctionOrNothing,
#             NumberOrNothing,
#             FunctionOrNothing,
#             FunctionOrNothing,
#             NumberOrNothing
#             } where {T<:Number, U<:Number}
#     }


@doc """
    ab_generalized(V, p)

The generalized form of α and β.

```
function ab_generalized(V, p)
    A, B, C, D, F, H = p
    (A + B * V) / (C + H * exp((V + D) / F))
end

```

"""
function ab_generalized(V, p)
    A, B, C, D, F, H = p
    (A + B * V) / (C + H * exp((V + D) / F))
end

@doc raw"""
    struct HHChannel{T <: Number, U <: Number}
        id::String
        g::T
        Vsteady::T
        alpha_m::Function
        beta_m::Function
        p::U
        alpha_h::Function
        beta_h::Function
        q::U
    end

General struct for a channel in the Hodgkin-Huxley model in the form:

```math
I_x = \bar{g}_xm^p_xh^q_x(V_mV_x)
```

The gating variables ``m`` and ``h`` are are given by `α_m`, `β_m`, `α_h` and `β_h` respectively in the following wa

```math
\begin{aligned}
\frac{dm}{dt} = \alpha(V_M)(1-m) - \beta(V_M)m \\
\frac{dh}{dt} = \alpha(V_M)(1-h) - \beta(V_M)h
\end{aligned}
```

"""
struct HHChannel{T <: Number, U <: Number}
    id::String
    g::T  # Max conductance
    Vsteady::T # Ion specific reversal potential
    alpha_m::Function # Alpha for gating variable m
    beta_m::Function # Beta for gating variable m
    p::U # Exponent of gating variable m
    alpha_h::Function # Alpha for gating variable h
    beta_h::Function # Beta for garing variable h
    q::U # Exponent for gating variable h
end

@doc """
    hhchannel(
    id,
    g::T,
    Vsteady::T;
    alpha_m=() -> oneunit(T),
    beta_m=() -> oneunit(T),
    p=zero(T),
    alpha_h=() -> oneunit(T),
    beta_h=() -> oneunit(T),
    q=zero(T)
    ) where T<:Number


Outer constructor of [`LaskaML.hh.HHChannel`](@ref).
If not specified, exponents of ``m`` and ``h`` will be 0 while unspecified α and β functions will be initialized as `() -> oneunit(T)`.

"""
function hhchannel(
    id,
    g::T,
    Vsteady::T;
    alpha_m=() -> oneunit(T),
    beta_m=() -> oneunit(T),
    p=zero(T),
    alpha_h=() -> oneunit(T),
    beta_h=() -> oneunit(T),
    q=zero(T)
    ) where T<:Number
    HHChannel(id,
              g,
              Vsteady,
              alpha_m,
              beta_m,
              p,
              alpha_h,
              beta_h,
              q)
end

# Interface for HHChannel

@doc """
    id(channel::HHChannel)

Get the identifier of a `HHChannel`.

"""
id(channel::HHChannel) = channel.id

@doc raw"""
    conductance(channel::HHChannel)

Get the maximum conductance, ``\bar{g}``, of a HHChannel.

"""
conductance(channel::HHChannel) = channel.g

@doc """
    alpha_m(channel::HHChannel) = channel.α_m

Get the function for ``\alpha_m`` bound to a HHChannel.

"""
alpha_m(channel::HHChannel) = channel.alpha_m

@doc """
    beta_m(channel::HHChannel)

Get the function for ``\beta_m`` bound to a HHChannel.

"""
beta_m(channel::HHChannel) = channel.beta_m

@doc """
    m_exponent(channel::HHChannel)

Get the exponent of ``m`` in a HHChannel

"""
m_exponent(channel::HHChannel) = channel.p

@doc """
    alpha_h(channel::HHChannel)

Get the function for ``\alpha_h`` bound to a HHChannel.

"""
alpha_h(channel::HHChannel) = channel.alpha_h

@doc """
    beta_h(channel::HHChannel)

Get the function for ``\beta_h`` bound to a HHChannel.

"""
beta_h(channel::HHChannel) = channel.beta_h

@doc """
    h_exponent(channel::HHChannel)

Get the exponent of ``h`` in a HHChannel.

"""
h_exponent(channel::HHChannel) = channel.q


# Parse a single channel

@doc raw"""
    function buildterm(variable::String, id::String, exponent::T) where T<:Number

Build a `String` representation of ``m^p`` or ``h^q``. If the exponent of either term is 0 it is not included.

"""
function buildterm(variable::String, id::String, exponent::T) where T<:Number
	iszero(exponent) ? "" : "($(variable)_$(id)^$(exponent))*"
end

@doc """
    parsechannel(channel::HHChannel)

Parse a `HHChannel` struct into an expression in the form of a `String` suitable for use in DifferentialEquations.

# Examples
```julia-repl
julia> chanK = LaskaML.hh.hhchannel(
               "K",
               35.0,
               -77.0,
               alpha_m=alpha_n,
               beta_m=beta_n,
               p=4,
               alpha_h=alpha_h,
               beta_h=beta_h,
               q=0
           )
LaskaML.hh.HHChannel{Float64, Int64}("K", 35.0, -77.0, alpha_n, beta_n, 4, alpha_h, beta_h, 0)

julia> LaskaML.hh.parsechannel(chanK)
"(gK*(m_K^4)*(v-E_K))"
```

"""
function parsechannel(channel::HHChannel)
    i = hh.id(channel)
    m_exp = buildterm("m", i, m_exponent(channel))
    h_exp = buildterm("h", i, h_exponent(channel))
    "(g_$(i)*$(m_exp)$(h_exp)(v-E_$(id(channel))))"
end

# Struct for holding an entire model

@doc """
    struct HHModel{T<:Number, U<:Number}
        channels::Vector{HHChannel{T, U}}

        function HHModel(channels::HHChannel{T,U}...) where {T<:Number, U<:Number}
	           new{T, U}([channels...])
        end
    end

Struct for holding several `HHChannels`, representing an (almost complete) model.

The channels are indexable using `model[i]`.

The struct can be initialized by calling its inner constructor and passing it the channels to be part of the model.
For convenience, the struct may also be initialized by "adding" `HHChannels` together using the `+` operator. The `+` operator may also be used to add `HHChannel`s to an existing `HHModel`

"""
struct HHModel{T<:Number, U<:Number}
    channels::Vector{HHChannel{T, U}}

    function HHModel(channels::HHChannel{T,U}...) where {T<:Number, U<:Number}
	      new{T, U}([channels...])
    end
end


channels(model::HHModel) = model.channels

# Interface functions for HHModel

# Basics

Base.length(X::HHModel) = length(channels(X))

# Indexing
function Base.getindex(X::HHModel, i::Int)
	  1 <= i <= length(X) || throw(BoundsError(X, i))
    return channels(X)[i]
end
function Base.setindex!(X::HHModel, i::Int, v::HHChannel)
	  channels(X)[i] = v
end
Base.firstindex(X::HHModel) = 1
Base.lastindex(X::HHModel) = length(X)

# Adding
Base.:+(x::HHModel, y::HHChannel...) = for c in y push!(channels(x), c) end
function Base.:+(x::HHChannel...)
    HHModel(x...)
end


end

