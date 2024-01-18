
module hh

const FunctionOrNothing = Union{Function, Nothing}
const NumberOrNothing = Union{Number, Nothing}

"""
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
    struct HHChannel{T <: FunctionOrNothing,U <: NumberOrNothing, V <: Number}
        g::V
        Vsteady::V
        α_m::T
        β_m::T
        p::U
        α_h::T
        β_h::T
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
struct HHChannel{T<:Number, U<:Number}
    identifier::String
    g::U  # Max conductance
    Vsteady::Number # Ion specific reversal potential
    α_m::FunctionOrNothing # Alpha for gating variable m
    β_m::FunctionOrNothing # Beta for gating variable m
    p::T # Exponent of gating variable m
    α_h::FunctionOrNothing # Alpha for gating variable h
    β_h::FunctionOrNothing # Beta for garing variable h
    q::T # Exponent for gating variable h
end

"""
    function hhchannel(g, α_m=nothing, β_m=nothing, p=nothing, α_h=nothing, β_h=nothing, q=nothing)

Outer constructor of the HHChannel. Mostly convenient for *leaky channels* and similar without gating variables.

"""
function hhchannel(identifier, g, Vsteady, α_m=nothing, β_m=nothing, p=nothing, α_h=nothing, β_h=nothing, q=nothing)
    HHChannel(identifier, g, Vsteady, α_m, β_m, p, α_h, β_h, q)
end

id(channel::HHChannel) = channel.identifier
conductance(channel::HHChannel) = channel.g
alpha_m(channel::HHChannel) = channel.α_m
beta_m(channel::HHChannel) = channel.β_m
m_exponent(channel::HHChannel) = channel.p
alpha_h(channel::HHChannel) = channel.α_h
beta_h(channel::HHChannel) = channel.β_h
h_exponent(channel::HHChannel) = channel.q


# Struct for holding an entire model

struct HHModel{T<:Number}
    channels::Vector{HHChannel}
    I::T
    C::T
end

channels(model::HHModel) = model.channels

function buildmodel(model::HHModel)
    parsed = String[]
	for (_, ch) in enumerate(channels(model))
    end
end

# Parse a single channel
function parsechannel(channel::HHChannel)
    "(g$(id(channel)) * (m$(id(channel))^$(m_exponent(channel))) * (h$(id(channel))^$(h_exponent(channel))) * (v-E$(id(channel))))"
end

end
