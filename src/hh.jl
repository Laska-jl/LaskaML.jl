
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
    struct HHChannel{G<:Number, V<:Number, P<:Number, Q<:Number}
        id::String
        g::G  # Max conductance
        Vrev::V # Ion specific reversal potential
        alpha_m::Union{Function, Nothing} # Alpha for gating variable m
        beta_m::Union{Function, Nothing} # Beta for gating variable m
        p::P # Exponent of gating variable m
        alpha_h::Union{Function, Nothing} # Alpha for gating variable h
        beta_h::Union{Function, Nothing} # Beta for garing variable h
        q::Q # Exponent for gating variable h
        m_inf::Union{Function, Nothing}
        h_inf::Union{Function, Nothing}
    end

General struct for a channel in the Hodgkin-Huxley model in the form:

```math
I_x = \bar{g}_xm^p_xh^q_x(V_m-V_x)
```

The gating variables ``m`` and ``h`` are are given by `α_m`, `β_m`, `α_h` and `β_h` respectively in the following wa

```math
\begin{aligned}
\frac{dm}{dt} = \alpha(V_M)(1-m) - \beta(V_M)m \\
\frac{dh}{dt} = \alpha(V_M)(1-h) - \beta(V_M)h
\end{aligned}
```

"""
struct HHChannel{G<:Number, V<:Number, P<:Number, Q<:Number}
    id::String
    g::G  # Max conductance
    Vrev::V # Ion specific reversal potential
    alpha_m::Union{Function, Nothing} # Alpha for gating variable m
    beta_m::Union{Function, Nothing} # Beta for gating variable m
    p::P # Exponent of gating variable m
    alpha_h::Union{Function, Nothing} # Alpha for gating variable h
    beta_h::Union{Function, Nothing} # Beta for garing variable h
    q::Q # Exponent for gating variable h
    m_inf::Union{Function, Nothing}
    h_inf::Union{Function, Nothing}
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


Outer constructor of [`HHChannel`](@ref).
If not specified, exponents of ``m`` and ``h`` will be 0 while unspecified α and β functions will be initialized as `nothing`.

"""
function hhchannel(
    id,
    g::G,
    Vrev::V;
    alpha_m=nothing,
    beta_m=nothing,
    p::P=zero(Int),
    alpha_h=nothing,
    beta_h=nothing,
    q::Q=zero(Int)
    ) where {G<:Number, V<:Number, P<:Number, Q<:Number}
    # m steady state
    if isnothing(alpha_m) && isnothing(beta_m)
        m_inf = nothing
    else
	    m_inf = (V) -> alpha_m(V) / (alpha_m(V) + beta_m(V))
    end
    # h steady state
    if isnothing(alpha_h) && isnothing(beta_h)
        h_inf = nothing
    else
	    h_inf = (V) -> alpha_h(V) / (alpha_h(V) + beta_h(V))
    end
    HHChannel(id,
              g,
              Vrev,
              alpha_m,
              beta_m,
              p,
              alpha_h,
              beta_h,
              q,
              m_inf,
              h_inf)
end

# Interface for HHChannel
Base.convert(::Type{HHChannel{G, V, P, Q}}, x::HHChannel{G, V, P, Q}) where {G<:Number, V<:Number, P<:Number, Q<:Number} = x

function Base.convert(::Type{HHChannel{G,V,P,Q}}, x::HHChannel) where {G<:Number, V<:Number, P<:Number, Q<:Number}
    HHChannel(
        id(x),
        G(conductance(x)),
        V(Vrev(x)),
        alpha_m(x),
        beta_m(x),
        P(m_exponent(x)),
        alpha_h,
        beta_h,
        Q(h_exponent(x))
    )
end
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
     Vrev(channel::HHChannel) = channel.Vsteady

Get the reversal potential.

"""
Vrev(channel::HHChannel) = channel.Vrev

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


"""
    m_inf(channel::HHChannel)

Get the function for steady state of ``m``.

"""
m_inf(channel::HHChannel) = channel.m_inf


"""
    h_inf(channel::HHChannel)

Get the function for steady state of ``h``.

"""
h_inf(channel::HHChannel) = channel.h_inf


# Parse a single channel

@doc raw"""
    function buildterm(variable::String, id::String, exponent::T) where T<:Number

Build a `String` representation of ``m^p`` or ``h^q``. If the exponent of either term is 0 it is not included.

"""
function buildterm(variable::String, id::String, exponent::T) where T<:Number
    if iszero(exponent)
        return ""
    elseif isone(exponent)
        return "$(variable)_$(id)*"
    else
        return "($(variable)_$(id)^$(exponent))*"
    end
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
    "-((g_$(i)*$(m_exp)$(h_exp)(V-E_$(id(channel)))))"
end

# Struct for holding an entire model

@doc """
    struct HHModel{G<:Number, V<:Number, P<:Number, Q<:Number}
        V_0::G
        I_0::G
        C_0::G
        channels::Vector{HHChannel{G,V,P,Q}}

        function HHModel(V_0, I_0, C_0, channels::HHChannel{G, V, P, Q}...) where {G<:Number, V<:Number, P<:Number, Q<:Number}
	        new{G,V,P,Q}(V_0, I_0, C_0, [channels...])
        end
    end

Struct for holding several `HHChannels`, representing an (almost complete) model.

The channels are indexable using `model[i]`.

The struct can be initialized by calling its inner constructor and passing it the channels to be part of the model.
For convenience, the struct may also be initialized by "adding" `HHChannels` together using the `+` operator. The `+` operator may also be used to add `HHChannel`s to an existing `HHModel`

"""
struct HHModel{G<:Number, V<:Number, P<:Number, Q<:Number}
    V_0::G
    I_0::G
    C_0::G
    channels::Vector{HHChannel{G,V,P,Q}}

    function HHModel(V_0, I_0, C_0, channels::HHChannel{G, V, P, Q}...) where {G<:Number, V<:Number, P<:Number, Q<:Number}
	      new{G,V,P,Q}(V_0, I_0, C_0, [channels...])
    end
end


channels(model::HHModel) = model.channels
V0(model::HHModel) = model.V_0
I0(model::HHModel) = model.I_0
C0(model::HHModel) = model.C_0

# Interface functions for HHModel

# Show

function Base.show(io::IO, ::MIME"text/plain", m::HHModel)
	print(
        io,
        "V_0: $(V0(m))\nI_O: $(I0(m))\nC_0: $(C0(m))\n$(join([repr(ch) for ch in channels(m)], "\n"))"
    )
end

function Base.show(io::IO, m::HHModel{G,V,P,Q}) where {G<:Number, V<:Number, P<:Number, Q<:Number}
	print(
        io,
        join([repr(ch) for ch in channels(m)], "\n")
    )
end

# Basics

@inline Base.length(X::HHModel) = length(channels(X))

# Indexing
@inline function Base.getindex(X::HHModel, i::Int)
	  1 <= i <= length(X) || throw(BoundsError(X, i))
    return channels(X)[i]
end
@inline function Base.setindex!(X::HHModel, i::Int, v::HHChannel)
	  channels(X)[i] = v
end
@inline Base.firstindex(X::HHModel) = 1
@inline Base.lastindex(X::HHModel) = length(X)

# Adding
@inline Base.:+(x::HHModel, y::HHChannel...) = for c in y push!(channels(x), c) end


@doc raw"""
    buildexpression(model::HHModel{G,V,P,Q}) where {G<:Number, V<:Number, P<:Number, Q<:Number}

Build a string representation of the expression:
```math
I = \bar{g}_xm_x^ph_x^q(V_M-V_x) + \bar{g}_ym_y^ph_y^q(V_M-V_y) + \bar{g}_zm_z^ph_z^q(V_M-V_z)...
```

Using all [`HHChannel`](@ref)s in the `model`. If the exponents ``p`` or ``q`` of a `HHChannel` is zero their respective term (``m`` or ``h``) will not be included.

# Examples

```julia-repl
julia> # Example potassium channel
           chanK = LaskaML.hh.hhchannel(
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
LaskaML.hh.HHChannel{Float64, Float64, Int64, Int64}("K", 35.0, -77.0, alpha_n, beta_n, 4, alpha_h, beta_h, 0)

julia> # Example sodium channel
           chanNa = LaskaML.hh.hhchannel(
               "Na",
               40.0,
               55.0,
               alpha_m=alpha_m,
               beta_m=beta_m,
               p=3,
               alpha_h=alpha_h,
               beta_h=beta_h,
               q=1
           )
LaskaML.hh.HHChannel{Float64, Float64, Int64, Int64}("Na", 40.0, 55.0, alpha_m, beta_m, 3, alpha_h, beta_h, 1)

julia> chanL = LaskaML.hh.hhchannel(
               "L",
               0.3,
               -65.0
           )
LaskaML.hh.HHChannel{Float64, Float64, Int64, Int64}("L", 0.3, -65.0, nothing, nothing, 0, nothing, nothing, 0)


julia> model = chanNa + chanK + chanL # Combine channels into a HHModel
LaskaML.hh.HHModel{Float64, Float64, Int64, Int64}(LaskaML.hh.HHChannel{Float64, Float64, Int64, Int64}[LaskaML.hh.HHChannel{Float64, Float64, Int64, Int64}("Na", 40.0, 55.0, alpha_m, beta_m, 3, alpha_h, beta_h, 1), LaskaML.hh.HHChannel{Float64, Float64, Int64, Int64}("K", 35.0, -77.0, alpha_n, beta_n, 4, alpha_h, beta_h, 0), LaskaML.hh.HHChannel{Float64, Float64, Int64, Int64}("L", 0.3, -65.0, nothing, nothing, 0, nothing, nothing, 0)])
nel{Float64, Float64, Int64, Int64}("L", 0.3, -65.0, nothing, nothing, 0, nothing, nothing, 0)

julia> LaskaML.hh.buildexpression(model) # Build String expression from model.
"(g_Na*(m_Na^3)*h_Na*(v-E_Na))+(g_K*(m_K^4)*(v-E_K))+(g_L*(v-E_L))"
```

"""
function buildexpression(model::HHModel{G,V,P,Q}) where {G<:Number, V<:Number, P<:Number, Q<:Number}
    res = Vector{String}(undef, length(model))
    for (n,ch) in enumerate(channels(model))
	      res[n] = parsechannel(ch)
    end
    "(" * join(res) * " + I)" * " / C"
end

@doc raw"""
     function dmdv(channel::HHChannel)

Returns a function giving dm/dv according to:

```math
\frac{dm}{dv} = \alpha(v)(1.0 - m) - \beta(v)m
```

"""
function dmdv(channel::HHChannel)
    am = alpha_m(channel)
    bm = beta_m(channel)
    function (in::Tuple{T, T}) where {T}
        V, m = in
        (am(V) * (oneunit(T) - m)) - (bm(V) * m)
    end
end

@doc raw"""
     function dhdv(channel::HHChannel)

Returns a function giving dh/dv according to:

```math
\frac{dh}{dv} = \alpha(v)(1.0 - h) - \beta(v)h
```

"""
function dhdv(channel::HHChannel)
    ah = alpha_h(channel)
    bh = beta_h(channel)
    function (in::Tuple{T, T}) where {T}
        V, h = in
        (ah(V) * (oneunit(T) - h)) - (bh(V) * h)
    end
end

"""
    buildd0(model::HHModel{G,V,P,Q}) where {G<:Number, V<:Number, P<:Number, Q<:Number}

Construct the initial conditions as well as the functions for ``frac{dm}{dv}`` and/or ``frac{dh}{dv}`` for each channel.

**Returns** a `Tuple` consisting of:

- A `Vector` of the functions mentioned above.
- A `Vector` of the initial conditions.
- A `Vector` containing `String`s identifying the functions.
"""
function buildd0(model::HHModel{G,V,P,Q}) where {G<:Number, V<:Number, P<:Number, Q<:Number}
    u0 = Vector{promote_type(G,V,P,Q)}(undef, 0)
    d = Vector{Function}(undef, 0)
    v = V0(model)
    push!(u0, v)
    textrepr = ["V"]
    for ch in channels(model)
        minf = m_inf(ch)
        if !isnothing(minf)
            push!(u0, minf(v))
            push!(textrepr, "m_" * id(ch))
            push!(d, dmdv(ch))
        end
        hinf = h_inf(ch)
        if !isnothing(hinf)
            push!(u0, hinf(v))
            push!(textrepr, "h_" * id(ch))
            push!(d, dhdv(ch))
        end
    end
    return d, u0, textrepr
end

"""
    buildparams(model::HHModel{G,V,P,Q}) where {G<:Number, V<:Number, P<:Number, Q<:Number}

returns the initial parameters of the `model` as well as a `Vector{String}` of their identifiers.
"""
function buildparams(model::HHModel{G,V,P,Q}) where {G<:Number, V<:Number, P<:Number, Q<:Number}
    outg = Vector{G}(undef, length(model))
    oute = Vector{G}(undef, length(model))
    outstrg = Vector{String}(undef, length(model))
    outstre = Vector{String}(undef, length(model))
    for (i, ch) in enumerate(channels(model))
        outg[i] = conductance(ch)
        outstrg[i] = "g_" * id(ch)
        oute[i] = Vrev(ch)
        outstre[i] = "E_" * id(ch)
    end
    return vcat(outg, oute, [C0(model), I0(model)]), vcat(outstrg, outstre, "C", "I")
end


"""
    buildmodel(model::HHModel{G,V,P,Q}) where {G<:Number, V<:Number, P<:Number, Q<:Number}

Returns a `Tuple` containing all parts necessary for constructing an `ODEProblem` of the `model` except the time span:

- A mutating function accepting the parameters `du, u, p, t`
- The initial conditions of the model
- The parameters of the model
"""
function buildmodel(model::HHModel{G,V,P,Q}) where {G<:Number, V<:Number, P<:Number, Q<:Number}
	par, ptxt = buildparams(model)
    funcs, d0, d0txt = buildd0(model)
    ps = Meta.parse(join(ptxt, ", "))
    d0s = Meta.parse(join(d0txt, ", "))

    exprs = Vector{Expr}(undef, length(funcs)+3)
    exprs[1] = :($ps = p)
    exprs[2] = :($d0s = u)
    exprs[3] = :(du[1] = $(Meta.parse(buildexpression(model))))
    for n in 1:length(funcs)
        ind = n + 3
        exprs[ind] = :(du[$(n+1)] = $(funcs[n])($:($(Meta.parse(d0txt[1])), $(Meta.parse(d0txt[n+1])))))
    end
    out = :()
    out.head = :function
    out.args = [:((du, u, p, t)), Expr(:block, exprs...)]
    return eval(out), d0, par
end


end
