module FastQuantities

using Unitful: Unitful, Dimension, Dimensions
import Base: ==, ≈, one, *, /, ^, inv, sqrt, log, convert, show, values

export AbstractDimensions, SIDimensions, 
       AbstractUnitfulScalar, UnitfulScalar,
       NoDims, 𝐓, 𝐋, 𝐌, 𝐈, 𝚯, 𝐍, 𝐉,
       dimexps, value, dimensions,
       @u_str

allmap(f, xs...) = allequal(length(x) for x in xs) && all(Iterators.map(f, xs...))

"""
    AbstractDimensions

Abstract type representing the dimensions of a physical quantity.
"""
abstract type AbstractDimensions <: Number end

"""
    dimexps(x::AbstractDimensions)

Represent `x` as a product of base dimensions raised to some powers and return a tuple of these powers.

For [`SIDimensions`](@ref), the base dimensions are
time, length, mass, current, temperature, amount, luminous intensity.
"""
dimexps(x::AbstractDimensions) = x.exps

show(io::IO, x::AbstractDimensions) = show(io, convert(Dimensions, x))

==(x::T, y::T) where T <: AbstractDimensions = dimexps(x) == dimexps(y)
≈(x::T, y::T) where T <: AbstractDimensions = allmap((ξ, η) -> ≈(ξ, η; atol = 2*sqrt(eps(Float32))),
                                                     dimexps(x), dimexps(y))
*(x::T, y::T) where T <: AbstractDimensions = T(dimexps(x) .+ dimexps(y))
/(x::T, y::T) where T <: AbstractDimensions = T(dimexps(x) .- dimexps(y))
^(x::T, y::Number) where T <: AbstractDimensions = T(dimexps(x) .* y)
^(x::T, y::Integer) where T <: AbstractDimensions = T(dimexps(x) .* y) # to resolve method ambiguity
^(x::T, y::Rational) where T <: AbstractDimensions = T(dimexps(x) .* y) # to resolve method ambiguity
inv(x::T) where T <: AbstractDimensions = T((-).(dimexps(x)))
sqrt(x::AbstractDimensions) = x^(1//2)

function log(x::AbstractDimensions)
    if x == one(x)
        return one(x)
    else
        throw(ArgumentError("log of dimensionful quantities is not supported
                until logarithmic quantities are implemented"))
    end
end

"""
    SIDimensions

Concrete type storing the dimensions of a physical quantity as 7 `Float32` exponents
of SI base units.
"""
struct SIDimensions <: AbstractDimensions
    exps::NTuple{7, Float32}
end

one(::Type{SIDimensions}) = SIDimensions(ntuple(_ -> 0, 7))

"""
    NoDims
    
Physical dimensions of a dimensionless quantity.
"""
const NoDims = one(SIDimensions)

const UnitfulBaseDims = (Unitful.𝐓, Unitful.𝐋, Unitful.𝐌, Unitful.𝐈, Unitful.𝚯, Unitful.𝐍, Unitful.𝐉)
const UnitfulBaseNames = UnitfulBaseDims .|> (((::Dimensions{D}) where D) -> D[1]) .|> Unitful.name

SIDimensions(x::Dimension{T}) where T = SIDimensions(Unitful.power(x) .* (Unitful.name(x) .== UnitfulBaseNames))
SIDimensions(x::Dimensions{D}) where D = prod(SIDimensions.(D), init=NoDims)

convert(::Type{Dimensions}, x::SIDimensions) = prod(UnitfulBaseDims .^ dimexps(x))
convert(::Type{SIDimensions}, x::Dimensions) = SIDimensions(x)

𝐓, 𝐋, 𝐌, 𝐈, 𝚯, 𝐍, 𝐉 = SIDimensions.(UnitfulBaseDims)

"""
    AbstractUnitfulScalar{TV, TD<:AbstractDimensions} <: Number

Abstract type representing a scalar physical quantity
with a numerical value of type `TV` and physical dimensions of type `TD`.
"""
abstract type AbstractUnitfulScalar{TV, TD<:AbstractDimensions} <: Number end

"""
    value(x::AbstractUnitfulScalar)

Get the numerical values of an [`AbstractUnitfulScalar`](@ref) in the default unit system (SI).

See also: [`values`](@ref values(::AbstractUnitfulScalar)),
[`dimensions`](@ref dimensions(::AbstractUnitfulScalar)).
"""
value(x::AbstractUnitfulScalar) = x.value

"""
    dimensions(x::AbstractUnitfulScalar)

Get the physical dimensions of an [`AbstractUnitfulScalar`](@ref).

See also: [`value`](@ref).
"""
dimensions(x::AbstractUnitfulScalar) = x.dims

"""
    values(x::AbstractUnitfulScalar)

Same as [`value`](@ref), for consistency with `AbstractUnitfulTensor`.

See also: [`dimensions`](@ref  dimensions(::AbstractUnitfulScalar)).
"""
values(x::AbstractUnitfulScalar) = value(x)

value(x::Number) = x
dimensions(x::Number) = NoDims
# values(x::Number) works as intended because of the values(itr) = itr definition in Base

show(io::IO, x::AbstractUnitfulScalar) = show(io, value(x) * Unitful.upreferred(convert(Dimensions, dimensions(x))))

"""
    UnitfulScalar{TV, TD<:AbstractDimensions} <: AbstractUnitfulScalar{TV, TD}

Concrete type representing a scalar physical quantity
with a numerical value of type `TV` and physical dimensions of type `TD`.
"""
struct UnitfulScalar{TV, TD<:AbstractDimensions} <: AbstractUnitfulScalar{TV, TD}
    value::TV
    dims::TD
end

macro u_str(str)
    quote
        unit = Unitful.@u_str($str)
        val = Unitful.ustrip(Unitful.upreferred(unit), 1. * unit)
        dims = convert(SIDimensions, Unitful.dimension(unit))
        UnitfulScalar(val, dims)
    end
end

end # of module
