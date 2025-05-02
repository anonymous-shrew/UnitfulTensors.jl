module FastQuantities

using Unitful: Unitful, Dimension, Dimensions
using DynamicQuantities: Quantity
import DynamicQuantities
import Base: ≈, log, convert, show, values, length, iterate, ndims, broadcastable

export AbstractDimensions, SIDimensions,
       AbstractUnitfulScalar, UnitfulScalar, Quantity,
       NoDims, 𝐓, 𝐋, 𝐌, 𝐈, 𝚯, 𝐍, 𝐉,
       dimexps, value, dimensions,
       @u_str

allmap(f, xs...) = allequal(length(x) for x in xs) && all(Iterators.map(f, xs...))

"""
    AbstractDimensions{R}

Abstract type representing the dimensions of a physical quantity. R is the type of dimensional exponents.

Alias for `DynamicQuantities.AbstractDimensions`.
"""
const AbstractDimensions = DynamicQuantities.AbstractDimensions

"""
    dimexps(x::AbstractDimensions)

Represent `x` as a product of base dimensions raised to some powers and return a tuple of these powers.

For [`SIDimensions`](@ref), the base dimensions are
time, length, mass, current, temperature, amount, luminous intensity.
"""
dimexps(x::AbstractDimensions) = x.exps

show(io::IO, x::AbstractDimensions{<:Number}) = show(io, convert(Dimensions, x)) # without <:Number this would overwrite a method from DynamicQuantities and break precompilation

≈(x::T, y::T) where T <: AbstractDimensions = allmap((ξ, η) -> ≈(ξ, η; atol = 2*sqrt(eps(Float32))),
                                                     dimexps(x), dimexps(y))
function log(x::AbstractDimensions)
    if x == one(x)
        return one(x)
    else
        throw(ArgumentError("log of dimensionful quantities is not supported
                until logarithmic quantities are implemented"))
    end
end

length(x::AbstractDimensions) = 1
iterate(x::AbstractDimensions) = (x, nothing)
iterate(x::AbstractDimensions, ::Any) = nothing
ndims(x::AbstractDimensions) = 0
broadcastable(x::AbstractDimensions) = Ref(x)

"""
    SIDimensions <: AbstractDimensions{Float32}

Concrete type storing the dimensions of a physical quantity as 7 `Float32` exponents
of SI base units.

Alias for `DynamicQuantities.Dimensions{Float32}`.
"""
const SIDimensions = DynamicQuantities.Dimensions{Float32}

const DynamicQuantitiesBaseNames = (:time, :length, :mass, :current, :temperature, :amount, :luminosity)
SIDimensions(exps::NTuple{7}) = SIDimensions(; Pair.(DynamicQuantitiesBaseNames, exps)...)
dimexps(x::SIDimensions) = getproperty.(x, DynamicQuantitiesBaseNames)

"""
    NoDims
    
Physical dimensions of a dimensionless quantity.
"""
const NoDims = one(SIDimensions)

const UnitfulBaseDims = (Unitful.𝐓, Unitful.𝐋, Unitful.𝐌, Unitful.𝐈, Unitful.𝚯, Unitful.𝐍, Unitful.𝐉)

SIDimensions(x::Dimensions) = convert(SIDimensions, x)

convert(::Type{Dimensions}, x::SIDimensions) = prod(UnitfulBaseDims .^ dimexps(x))

𝐓, 𝐋, 𝐌, 𝐈, 𝚯, 𝐍, 𝐉 = SIDimensions.(UnitfulBaseDims)

"""
    AbstractUnitfulScalar{TV, TD<:AbstractDimensions} <: Number

Abstract type representing a scalar physical quantity
with a numerical value of type `TV` and physical dimensions of type `TD`.

Alias for `DynamicQuantities.AbstractQuantity{TV, TD} where {TV, TD<:AbstractDimensions}`.
"""
const AbstractUnitfulScalar = DynamicQuantities.AbstractQuantity{TV, TD} where {TV, TD<:AbstractDimensions}

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
dimensions(x::AbstractUnitfulScalar) = x.dimensions

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
show(io::IO, x::AbstractUnitfulScalar{<:Real}) = show(io, value(x) * Unitful.upreferred(convert(Dimensions, dimensions(x)))) # resolving ambiguity with DynamicQuantities

"""
    UnitfulScalar{TV<:Number, TD<:AbstractDimensions} <: AbstractUnitfulScalar{TV, TD}

Concrete type representing a scalar physical quantity
with a numerical value of type `TV` and physical dimensions of type `TD`.

Alias for `DynamicQuantities.Quantity`.
"""
const UnitfulScalar = DynamicQuantities.Quantity

UnitfulScalar{T, D}(val, dims) where {T<:Number,D<:AbstractDimensions} = UnitfulScalar(convert(T, val), convert(D, dims))

macro u_str(str)
    quote
        unit = Unitful.@u_str($str)
        val = Unitful.ustrip(Unitful.upreferred(unit), 1. * unit)
        dims = convert(SIDimensions, Unitful.dimension(unit))
        UnitfulScalar(val, dims)
    end
end

# fixes the printing of type aliases, https://github.com/JuliaLang/julia/issues/40448
Base.modulesof!(s::Set{Module}, x::Type{<:Union{SIDimensions, AbstractUnitfulScalar}}) = (push!(s, @__MODULE__); s)

end # of module
