
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

"""
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

General struct for a channel in the Hodgkin-Huxley model.

The general expression of the channel is:

```math
I_x = \bar{g}_xm^p_xh^q_x(V_mV_x)
```

The gating variables ``m`` and ``h`` are are given by `α_m`, `β_m`, `α_h` and `β_h` respectively in the following way:

```math
\begin{aligned}
\frac{dm}{dt} = \alpha(V_M)(1-m) - \beta(V_M)m \\
\frac{dh}{dt} = \alpha(V_M)(1-h) - \beta(V_M)h
\end{aligned}
```

"""
struct HHChannel{T <: FunctionOrNothing,U <: NumberOrNothing, V <: Number}
    g::V  # Max conductance
    Vsteady::V # Steady-state voltage
    α_m::T # Alpha for gating variable m
    β_m::T # Beta for gating variable m
    p::U # Exponent of gating variable m
    α_h::T # Alpha for gating variable n
    β_h::T # Beta for garing variable n
    q::U # Exponent for gating variable h
end

"""
    function hhchannel(g, α_m=nothing, β_m=nothing, p=nothing, α_h=nothing, β_h=nothing, q=nothing)

Outer constructor of the HHChannel. Mostly convenient for *leaky channels* and similar without gating variables.

"""
function hhchannel(g, Vsteady, α_m=nothing, β_m=nothing, p=nothing, α_h=nothing, β_h=nothing, q=nothing)
    HHChannel(g, Vsteady, α_m, β_m, p, α_h, β_h, q)
end

end
