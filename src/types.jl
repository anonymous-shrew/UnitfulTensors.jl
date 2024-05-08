"""
    AbstractAxisDimensions{T<:AbstractDimensions} <: AbstractVector{T}

Represents the physical dimensions of an [`AbstractUnitfulTensor`](@ref) along some axis.

Used as a building block of [`AbstractAxesDimensions`](@ref).
"""
abstract type AbstractAxisDimensions{T<:AbstractDimensions} <: AbstractVector{T} end

"""
    dimsvec(a::AbstractAxisDimensions)

Convert [`AbstractAxisDimensions`](@ref) to a `Vector` of [`AbstractDimensions`](@ref).

For [`AxisDimensions`](@ref), `dimsvec` returns the `normdims` field.

!!! warning
    Don't mutate it, or it will affect any [`UnitfulTensor`](@ref)
    that references this [`AxisDimensions`](@ref) internally. Instead, create a copy if necessary.
"""
dimsvec(a::AbstractAxisDimensions) = a.normdims

"""
    match(a::AbstractAxisDimensions, b::AbstractAxisDimensions)

Check if two [`AbstractAxisDimensions`](@ref) can be multiplied.

This function is used for dimensions checking in [`AbstractUnitfulMatrix`](@ref) multiplication.
"""
match(a::AbstractAxisDimensions, b::AbstractAxisDimensions) = allmap(match, a, b)
match(a::AbstractDimensions, b::AbstractDimensions) = inv(a) == b

"""
    AxisDimensions{T<:AbstractDimensions} <: AbstractAxisDimensions{T}

Represents the physical dimensions of a [`UnitfulTensor`](@ref) along some axis.

Used as a building block of [`AxesDimensions`](@ref).
"""
struct AxisDimensions{T<:AbstractDimensions} <: AbstractAxisDimensions{T}
    normdims::Vector{T}
    
    AxisDimensions{T}(dims) where T = new{T}(isempty(dims) ? dims : (dims ./ first(dims)))
    global _unsafe_AxisDimensions(dims::AbstractVector{T}) where T = new{T}(dims)
end

"""
    AxisDimensions(dims::AbstractVector; normalize = true)
    AxisDimensions(dims::AbstractVector{<:AbstractDimensions}; normalize = true)

Convert an `AbstractVector` of [`AbstractDimensions`](@ref) to [`AxisDimensions`](@ref).

By default, the dimensions are scaled by a scalar factor
so that the first element of [`AxisDimensions`](@ref) is [`NoDims`](@ref).
If you are absolutely sure that dims already satisfy this condition,
you may set `normalize = false` for performance.
"""
function AxisDimensions(dims::AbstractVector{T}; normalize = true) where T <: AbstractDimensions
    return normalize ? AxisDimensions{T}(dims) : _unsafe_AxisDimensions(dims)
end

function AxisDimensions(dims::AbstractVector; normalize = true)
    SIdims = convert(Vector{SIDimensions}, dims)
    return AxisDimensions(SIdims; normalize = normalize)
end

AxisDimensions(dims::AxisDimensions; kwargs...) = dims

"""
    AbstractAxesDimensions{N, T<:AbstractDimensions} <: AbstractArray{T, N}

Represents the physical dimensions of an `N`-dimensional [`AbstractUnitfulTensor`](@ref).
"""
abstract type AbstractAxesDimensions{N, T<:AbstractDimensions} <: AbstractArray{T, N} end

"""
    AbstractVectorDimensions{T<:AbstractDimensions}

Alias for [`AbstractAxesDimensions{1, T}`](@ref).
"""
const AbstractVectorDimensions{T<:AbstractDimensions} = AbstractAxesDimensions{1, T}

"""
    AbstractMatrixDimensions{T<:AbstractDimensions}

Alias for [`AbstractAxesDimensions{2, T}`](@ref).
"""
const AbstractMatrixDimensions{T<:AbstractDimensions} = AbstractAxesDimensions{2, T}
const AbstractVecOrMatDimensions{T<:AbstractDimensions} = Union{AbstractVectorDimensions{T},
                                                                AbstractMatrixDimensions{T}}

const AbstractDimensionsOrAxesDimensions = Union{AbstractDimensions, AbstractAxesDimensions}

"""
    normdims(x::AbstractAxesDimensions)
    normdims(x::AbstractAxisDimensions)
    normdims(x::AbstractDimensions)
    normdims(x::AbstractUnitfulTensor)

Get a tuple of [`AbstractAxisDimensions`](@ref) representing the physical dimensions of `x` along each axis.

[`AbstractAxisDimensions`](@ref) are normalized so that their first element is [`NoDims`](@ref).
The scalar normalization factor can be retrieved via [`dimscale(x)`](@ref).

For example, if `x` `isa` [`AbstractMatrixDimensions`](@ref), `normdims` returns its row and column dimensions,
and the elements of `x` satisfy `x[i, j] == dimscale(x) * normdims(x)[1][i] * normdims(x)[2][j]`.

[`AbstractAxisDimensions`](@ref) and [`AbstractDimensions`](@ref) are treated as
1- and 0-dimensional [`AbstractAxesDimensions`](@ref), respectively.
"""
normdims(x::AbstractAxesDimensions) = x.normdims

"""
    dimscale(x::AbstractAxesDimensions)
    dimscale(x::AbstractAxisDimensions)
    dimscale(x::AbstractDimensions)
    dimscale(x::AbstractUnitfulTensor)

Get the first element of `x`.

[`AbstractAxesDimensions`](@ref) are factorized into an overall scalar factor, retrieved by `dimscale`,
and zero or more [`AbstractAxisDimensions`](@ref), which represent the physical dimensions along each axis.
[`AbstractAxisDimensions`](@ref) are normalized so that their first element is [`NoDims`](@ref)
and can be retrieved via [`normdims(x)`](@ref). Normalization ensures that the factorization is unique.

For example, if `x` `isa` [`AbstractMatrixDimensions`](@ref),
the elements of `x` satisfy `x[i, j] == dimscale(x) * normdims(x)[1][i] * normdims(x)[2][j]`.

[`AbstractAxisDimensions`](@ref) and [`AbstractDimensions`](@ref) are treated as
1- and 0-dimensional [`AbstractAxesDimensions`](@ref), respectively.
"""
dimscale(x::AbstractAxesDimensions) = x.scale

normdims(x::AbstractAxisDimensions) = (x, )
dimscale(x::AbstractAxisDimensions{T}) where T = one(T)

normdims(x::AbstractDimensions) = ()
dimscale(x::AbstractDimensions) = x

"""
    dimsplat(x)

Return [`dimscale(x)`](@ref) and all [`normdims(x)`](@ref) as a single tuple.
"""
dimsplat(x) = (dimscale(x), normdims(x)...)

==(A::T, B::T) where T <: AbstractAxesDimensions = ==(dimscale.((A, B))...) && ==(normdims.((A, B))...)

"""
    AxesDimensions{N, T<:AbstractDimensions} <: AbstractAxesDimensions{N, T}

Represents the physical dimensions of an `N`-dimensional [`UnitfulTensor`](@ref).
"""
struct AxesDimensions{N, T<:AbstractDimensions} <: AbstractAxesDimensions{N, T}
    normdims::NTuple{N, AxisDimensions{T}}
    scale::T
end

"""
    AxesDimensions(dims, scale; normalize = true)
    AxesDimensions(dims; normalize = true)

Assemble [`AxesDimensions`](@ref) from the physical dimensions along each axis (`dims`)
and an overall scalar factor (`scale`).

The inverse operation of factorizing [`AxesDimensions`](@ref) into `dims` and `scale`
can be performed via [`normdims`](@ref) and [`dimscale`](@ref).

`dims` can be a tuple of [`AxisDimensions`](@ref) or a tuple of `Vector`s of [`AbstractDimensions`](@ref).
If `scale == NoDims`, it can be omitted.

By default, `dims` are scaled by a scalar factor so that their first element is [`NoDims`](@ref).
If you are absolutely sure that dims already satisfy this condition,
you may set `normalize = false` for performance.
"""
function AxesDimensions(dims, scale; normalize = true)
    # Because AxisDimensions.(dims; normalize = normalize) is slow
    f(x) = AxisDimensions(x; normalize = normalize)
    return AxesDimensions(f.(dims), scale)
end

AxesDimensions(dims::Tuple; kwargs...) = AxesDimensions(dims, one(eltype(first(dims))); kwargs...)
AxesDimensions(::Tuple{}; kwargs...) = AxesDimensions((), NoDims; kwargs...)

"""
    AxesDimensions(A::AbstractArray{<:AbstractDimensions})

Convert an `AbstractArray` of physical dimensions to [`AxesDimensions`](@ref),
which enables efficient operations on the dimensions of unitful arrays.
"""
AxesDimensions(A::AbstractArray{<:AbstractDimensions}) = _get_single_index(A, eachindex(IndexCartesian(), A))
# _get_single_index is defined in indexing.jl

"""
    nodims(size...)

Create [`AxesDimensions`](@ref) representing the physical dimensions of a
dimensionless array of size `size`.
# Examples
```jldoctest
julia> nodims(2, 3)
2×3 AxesDimensions{2, SIDimensions}:
 NoDims  NoDims  NoDims
 NoDims  NoDims  NoDims

julia> UnitfulTensor([1 2 3; 4 5 6], 𝐍 * nodims(2, 3))
2×3 UnitfulMatrix{Int64, SIDimensions, Matrix{Int64}, AxesDimensions{2, SIDimensions}}:
 1 mol  2 mol  3 mol
 4 mol  5 mol  6 mol
```
"""
@inline nodims(size...) = AxesDimensions(fill.(NoDims, size), NoDims)

"""
    AbstractUnitfulTensor{T<:AbstractUnitfulScalar, N} <: AbstractArray{T, N}

An `N`-dimensional array of unitful quantities of type `T` that can be used in linear or multilinear algebra.

This requirement forces its physical dimensions,
represented by [`AbstractAxesDimensions`](@ref),
to factorize into a tensor product of dimensions along each axis,
represented by [`AbstractAxisDimensions`](@ref).
For example, the physical dimensions of an [`AbstractUnitfulMatrix`](@ref)
are the product of its row and column dimensions.
"""
abstract type AbstractUnitfulTensor{T<:AbstractUnitfulScalar, N} <: AbstractArray{T, N} end

"""
    AbstractUnitfulVector{T<:AbstractUnitfulScalar}

Alias for [`AbstractUnitfulTensor{T, 1}`](@ref).
"""
const AbstractUnitfulVector{T<:AbstractUnitfulScalar} = AbstractUnitfulTensor{T, 1}

"""
    AbstractUnitfulMatrix{T<:AbstractUnitfulScalar}

Alias for [`AbstractUnitfulTensor{T, 2}`](@ref).
"""
const AbstractUnitfulMatrix{T<:AbstractUnitfulScalar} = AbstractUnitfulTensor{T, 2}
const AbstractUnitfulVecOrMat{T<:AbstractUnitfulScalar} = Union{AbstractUnitfulVector{T},
                                                                AbstractUnitfulMatrix{T}}

const AbstractUnitfulScalarOrTensor = Union{AbstractUnitfulScalar, AbstractUnitfulTensor}

"""
    values(A::AbstractUnitfulTensor)

Get the numerical values of an [`AbstractUnitfulTensor`](@ref) in the default unit system (SI).

See also: [`dimensions`](@ref dimensions(::AbstractUnitfulTensor)).
"""
values(A::AbstractUnitfulTensor) = A.values
"""
    dimensions(A::AbstractUnitfulTensor)

Get the physical dimensions of an [`AbstractUnitfulTensor`](@ref).

See also: [`values`](@ref values(::AbstractUnitfulTensor)).
"""
dimensions(A::AbstractUnitfulTensor) = A.dims
dimensions(A::AbstractArray) = nodims(size(A)...)

normdims(x::AbstractUnitfulTensor) = normdims(dimensions(x))
dimscale(x::AbstractUnitfulTensor) = dimscale(dimensions(x))

"""
    UnitfulTensor{N, TV, TD<:AbstractDimensions,
                     V <: Union{AbstractArray{TV, N}, AbstractQ{TV}},
                     D <: AbstractAxesDimensions{N, TD}} <:
                     AbstractUnitfulTensor{UnitfulScalar{TV, TD}, N}

An `N`-dimensional array of unitful quantities that can be used in linear or multilinear algebra.

This requirement forces its physical dimensions,
represented by [`AxesDimensions`](@ref),
to factorize into a tensor product of dimensions along each axis,
represented by [`AxisDimensions`](@ref).
For example, the physical dimensions of a [`UnitfulMatrix`](@ref)
are the product of its row and column dimensions.

`UnitfulTensor` stores its numerical values as an array of type `V`
and its dimensions as an [`AbstractAxesDimensions`](@ref) of type `D`.
They can be obtained using [`values`](@ref values(::AbstractUnitfulTensor))
and [`dimensions`](@ref dimensions(::AbstractUnitfulTensor)).
"""
struct UnitfulTensor{N, TV, TD <: AbstractDimensions,
                     V <: Union{AbstractArray{TV, N}, AbstractQ{TV}},
                     D <: AbstractAxesDimensions{N, TD}} <:
                     AbstractUnitfulTensor{UnitfulScalar{TV, TD}, N}    
    values::V
    dims::D
end

"""
    UnitfulVector{TV, TD<:AbstractDimensions}

Alias for [`UnitfulTensor{1, TV, TD}`](@ref).
"""
const UnitfulVector{TV, TD<:AbstractDimensions} = UnitfulTensor{1, TV, TD}

"""
    UnitfulMatrix{TV, TD<:AbstractDimensions}

Alias for [`UnitfulTensor{2, TV, TD}`](@ref).
"""
const UnitfulMatrix{TV, TD<:AbstractDimensions} = UnitfulTensor{2, TV, TD}
const UnitfulVecOrMat{TV, TD<:AbstractDimensions} = Union{UnitfulVector{TV, TD},
                                                          UnitfulMatrix{TV, TD}}
const UnitfulScalarOrTensor = Union{UnitfulScalar, UnitfulTensor}

"""
    UnitfulTensor(A::AbstractArray)

Convert an `AbstractArray` of [`UnitfulScalar`](@ref)s or `Number`s to a [`UnitfulTensor`](@ref),
which enables efficient operations on the physical dimensions.

The physical dimensions of `A` must factorize into a tensor product of dimensions along each axis
(e. g., row and column dimensions of a matrix). Other arrays cannot be used in linear or
multilinear algebra and will throw an error when attempting to convert them to a [`UnitfulTensor`](@ref).

# Examples:
```jldoctest
julia> A = UnitfulTensor([1.0     2.0u"s^-1"
                          3.0u"J" 4.0u"W"   ])
2×2 UnitfulMatrix{Float64, SIDimensions, Matrix{Float64}, AxesDimensions{2, SIDimensions}}:
             1.0         2.0 s^-1
 3.0 kg m^2 s^-2  4.0 kg m^2 s^-3

julia> V = values(A)
2×2 Matrix{Float64}:
 1.0  2.0
 3.0  4.0

julia> D = dimensions(A)
2×2 AxesDimensions{2, SIDimensions}:
     NoDims        𝐓^-1
 𝐋^2 𝐌 𝐓^-2  𝐋^2 𝐌 𝐓^-3

julia> A == UnitfulTensor(V, D)
true

julia> inv(A)
2×2 UnitfulMatrix{Float64, SIDimensions, Matrix{Float64}, AxesDimensions{2, SIDimensions}}:
  -2.0   1.0 s^2 kg^-1 m^-2
 1.5 s  -0.5 s^3 kg^-1 m^-2

julia> inv(V)
2×2 Matrix{Float64}:
 -2.0   1.0
  1.5  -0.5
```
"""
function UnitfulTensor(A::AbstractArray)
    vals = value.(A); dims = dimensions.(A)
    return UnitfulTensor(vals, AxesDimensions(dims))
end

"""
    units_off()
Switch from [`UnitfulScalar`](@ref)s and [`UnitfulTensor`](@ref)s to plain numbers and arrays of numbers,
as if this package was not loaded.

If you are concerned with the overhead of using unitful quantities:
- Run your code on a small-scale problem to check the units
- Call `units_off()` immediately after `using $(@__MODULE__)`
- Run your code on the actual large-scale problem

If you used [`UnitfulScalar(...)`](@ref UnitfulScalar), [`UnitfulTensor(...)`](@ref UnitfulTensor),
and the `number * u"unit"` syntax for constructing unitful quantities, `units_off()` should eliminate
all runtime overhead. Quantities are still converted to the default units (SI) upon construction,
but then the units are stripped and the rest of the execution proceeds with plain numbers.

However, inner constructors like `UnitfulScalar{Float64, SIDimensions}(1., 𝐓)` will still produce
unitful quantities, so don't use them.

There is no `units_on` switch at present. You will have to restart Julia if you need units again.

# Examples
```julia
using UnitfulTensors, BenchmarkTools

13u"cm" |> println
typeof(13u"cm") |> println
u = UnitfulTensor([1u"eV", 2u"ns", 3u"μF"])
v = [1.602176634e-19, 1e-9, 3e-6]
@btime \$v * \$v'
@btime \$u * \$u'

units_off()
println()

13u"cm" |> println
typeof(13u"cm") |> println
u = UnitfulTensor([1u"eV", 2u"ns", 3u"μF"])
v = [1.602176634e-19, 1e-9, 3e-6]
@btime \$v * \$v'
@btime \$u * \$u'

# output

0.13 m
UnitfulScalar{Float64, SIDimensions}
  73.610 ns (1 allocation: 128 bytes)
  80.519 ns (1 allocation: 128 bytes)

0.13
Float64
  73.561 ns (1 allocation: 128 bytes)
  73.077 ns (1 allocation: 128 bytes)

3×3 Matrix{Float64}:
 2.56697e-38  3.20435e-28  4.80653e-25
 3.20435e-28  4.0e-18      6.0e-15
 4.80653e-25  6.0e-15      9.0e-12
```
"""
function units_off()
    for T in (:UnitfulScalar, :UnitfulTensor)
        M = @__MODULE__
        @eval delete_method.(methods($T, [$M, $M.FastQuantities]))
        @eval $T(args...) = args[1]
    end
end
