module Infinities

import Base: angle, isone, iszero, isinf, isfinite, isnan, isreal, abs, one, oneunit, zero, isless, isequal, inv,
                +, -, *, /, ^, ==, <, ≤, >, ≥, fld, cld, div, mod, rem, divrem, min, max,
                sign, signbit, isapprox,
                string, show, promote_rule, convert, getindex, tryparse, conj,
                isinteger, round, floor, ceil, trunc, float,
                Bool, Integer

export ∞,  ℵ₀,  ℵ₁, RealInfinity, ComplexInfinity, InfiniteCardinal, NotANumber, PositiveInfinity, NegativeInfinity
# The following is commented out for now to avoid conflicts with Infinity.jl
# export Infinity

"""
    NotANumber()

represents something that is undefined, for example, `0 * ∞`.

Every float type has a `NaN` of its own. This one belongs to none of them.
"""
struct NotANumber <: Real end

(::Type{T})(::NotANumber) where {T<:AbstractFloat} = T(NaN)
float(::NotANumber) = NaN
Base.hash(::NotANumber, h::UInt)::UInt = hash(NaN, h)


"""
   Infinity()

represents the positive real infinite.
"""
struct Infinity <: Real end

const ∞ = Infinity()

show(io::IO, ::Infinity) = print(io, "∞")
string(::Infinity) = "∞"

_convert(::Type{Float64}, ::Infinity) = Inf64
_convert(::Type{Float32}, ::Infinity) = Inf32
_convert(::Type{Float16}, ::Infinity) = Inf16
_convert(::Type{T}, ::Infinity) where {T<:Real} = convert(T, Inf)::T
(::Type{T})(x::Infinity) where {T<:Real} = _convert(T, x)

sign(y::Infinity) = 1
angle(x::Infinity) = 0
signbit(::Infinity) = false

one(::Type{Infinity}) = 1
oneunit(::Type{Infinity}) = 1
oneunit(::Infinity) = 1
zero(::Infinity) = 0
zero(::Type{Infinity}) = 0

abstract type RealInfinity <: Real end
struct PositiveInfinity <: RealInfinity end
struct NegativeInfinity <: RealInfinity end

signbit(::PositiveInfinity) = false
signbit(::NegativeInfinity) = true
one(::RealInfinity) = 1.0

RealInfinity() = PositiveInfinity()
RealInfinity(::Infinity) = PositiveInfinity()
RealInfinity(x::RealInfinity) = x
RealInfinity(x::Bool) = ifelse(x, NegativeInfinity(), PositiveInfinity())
PositiveInfinity(::Infinity) = PositiveInfinity() # otherwise the generic `(::Type{T})(::Infinity) where T<:Real` would route through `Inf`

_convert(::Type{Float16}, x::RealInfinity) = sign(x)*Inf16
_convert(::Type{Float32}, x::RealInfinity) = sign(x)*Inf32
_convert(::Type{Float64}, x::RealInfinity) = sign(x)*Inf64
_convert(::Type{T}, x::RealInfinity) where {T<:Real} = sign(x)*convert(T, Inf)
(::Type{T})(x::RealInfinity) where {T<:Real} = _convert(T, x)

for Typ in (RealInfinity, Infinity)
    @eval Bool(x::$Typ) = throw(InexactError(:Bool, Bool, x)) # ambiguity fix
end

sign(y::RealInfinity) = 1-2signbit(y)
angle(x::RealInfinity) = π*signbit(x)

string(y::RealInfinity) = signbit(y) ? "-∞" : "+∞"
show(io::IO, y::RealInfinity) = print(io, string(y))

Base.to_index(i::RealInfinity) = convert(Integer, i)

one(::Type{RealInfinity}) = 1.0
oneunit(::Type{RealInfinity}) = 1.0
oneunit(::RealInfinity) = 1.0
zero(::RealInfinity) = 0.0
zero(::Type{RealInfinity}) = 0.0


#######
# ComplexInfinity
#######

"""
    ComplexInfinity(halfturns)

represents an infinity in the complex plane, pointing in the direction `π * halfturns`.

The direction wraps at a full turn `τ = 2π`, so `ComplexInfinity(-1/2)` and
`ComplexInfinity(3/2)` are one and the same value.

See also [`reinterpret`](@ref).
"""
struct ComplexInfinity <: Number
    turns::UInt64
    ComplexInfinity(halfturns::Real) = new(_turns(halfturns))
    global _direction(turns::UInt64) = new(turns)
end

# A full turn τ fills the range, so a half turn is 2^63 units and neither conversion rounds.
@inline _turns(halfturns::Real) = round(UInt64, mod(halfturns, 2) * 0x1p63)
_turns(halfturns::Rational) = round(UInt64, mod(halfturns, 2) * big(2)^63)
@inline _halfturns(x::ComplexInfinity) = x.turns / 0x1p63

"""
    reinterpret(UInt64, x::ComplexInfinity)
    reinterpret(ComplexInfinity, turns::UInt64)

convert between an infinity and its direction as a count of `2^-64` turns.

The count wraps at a full turn, so every `UInt64` names a direction and every direction has
one count. `0x0` points along the positive real axis and `0x8000000000000000` along the
negative one.
"""
Base.reinterpret(::Type{UInt64}, x::ComplexInfinity) = x.turns
Base.reinterpret(::Type{ComplexInfinity}, turns::UInt64) = _direction(turns)

const _HALFTURN = 0x8000000000000000 # the negative real axis, which is its own negation

ComplexInfinity() = ComplexInfinity(false)
ComplexInfinity(::Infinity) = ComplexInfinity()
ComplexInfinity(x::RealInfinity) = ComplexInfinity(signbit(x))
ComplexInfinity(x::ComplexInfinity) = x

signbit(y::ComplexInfinity) = y.turns == _HALFTURN
isreal(y::ComplexInfinity) = iszero(y.turns) || signbit(y)

# `Base` converts a `Complex` to a `Real` the same way, and throws the same error off the axis.
RealInfinity(x::ComplexInfinity) = isreal(x) ? RealInfinity(signbit(x)) :
                                   throw(InexactError(:RealInfinity, RealInfinity, x))

convert(::Type{ComplexInfinity}, ::Infinity) = ComplexInfinity()
convert(::Type{ComplexInfinity}, x::RealInfinity) = ComplexInfinity(x)


sign(y::ComplexInfinity) = cispi(_halfturns(y))
angle(x::ComplexInfinity) = _halfturns(x) * π
abs(::ComplexInfinity) = ∞
conj(y::ComplexInfinity) = _direction(-y.turns)

# An exact zero has to stay finite, `Inf * 0` being a `NaN`.
@inline _ray(c) = iszero(c) ? c : copysign(Inf, c)
# `Complex` reaches only the eight rays of its two saturating parts, so the direction lands on the nearest of them.
function float(x::ComplexInfinity)
    s, c = sincospi(_halfturns(x))
    complex(_ray(c), _ray(s))
end

show(io::IO, x::ComplexInfinity) = print(io, "exp($(_halfturns(x))*im*π)∞")

one(::Type{ComplexInfinity}) = one(ComplexF64)
oneunit(::Type{ComplexInfinity}) = oneunit(ComplexF64)
oneunit(::ComplexInfinity) = oneunit(ComplexF64)
zero(::ComplexInfinity) = zero(ComplexF64)
zero(::Type{ComplexInfinity}) = zero(ComplexF64)


# `isequal` implies equal hashes, so the infinities have to hash like the float
# infinities they compare equal to. The interface requires implementing `hash(x, h::UInt)`.

Base.hash(::Infinity, h::UInt)::UInt = hash(Inf, h)
Base.hash(::PositiveInfinity, h::UInt)::UInt = hash(Inf, h)
Base.hash(::NegativeInfinity, h::UInt)::UInt = hash(-Inf, h)

# The two real directions have to hash like the real infinities they compare equal to.
function Base.hash(x::ComplexInfinity, h::UInt)::UInt
    iszero(x.turns) && return hash(Inf, h)
    x.turns == _HALFTURN && return hash(-Inf, h)
    hash(ComplexInfinity, hash(x.turns, h))
end


include("cardinality.jl")
include("interface.jl")
include("compare.jl")
include("algebra.jl")
include("ambiguities.jl")
end # module
