# Guide

## Creating UnitfulScalars

Unitful scalars can be defined using the u-prefixed non-standard string literals,
e. g.

```julia
Ry = 1.21104020061 * u"H*Mg*N*K*Pa*S*Tb*Tl*W/(Ba*C*Es*P*Pb*Pm*Ra*V*Yb)"
ħ = 1054.571818 * u"P/PP*pP*ppm*m*M/MM*Mm*mm"
```

They are parsed by [`Unitful.jl`](https://painterqubits.github.io/Unitful.jl/stable/)
and then converted to a [`UnitfulScalar`](@ref).

Don't try to use units defined on interval or logarithmic scales,
such as °C or dBm. Use K and mW instead.

## Creating UnitfulTensors

There are two ways to create [`UnitfulTensor`](@ref)s.
You can first construct an array of [`UnitfulScalar`](@ref)s

```julia
[1u"s" 2u"s"
 3u"s" 4u"s"]
```

and then
convert it to a [`UnitfulTensor`](@ref)

```julia
UnitfulTensor([1u"s" 2u"s"
               3u"s" 4u"s"])
```

Or you can construct an array of numbers and [`AxesDimensions`](@ref) separately
and then combine them into a [`UnitfulTensor`](@ref):

```julia
vals = [1. 2.
	3. 4.]
dims = 𝐓 * nodims(2, 2)
UnitfulTensor(vals, dims)
```

In the latter case, numerical values should be given in SI units.

Unitless quantities and arrays thereof don't have to be wrapped in
a [`UnitfulScalar`](@ref) or [`UnitfulTensor`](@ref),
but trying to do mixed unitful/unitless arithmetic with something
other than `Array`s (e. g., sparse matrices or `StaticArray`s)
can result in method ambiguities. They are easily fixed if necessary.

## Creating AxesDimensions

[`AxesDimensions`](@ref) represent the physical dimensions of a
[`UnitfulTensor`](@ref) and can be constructed as a tensor product
of one-dimensional [`AxesDimensions`](@ref):

```julia
AxesDimensions([𝐓, 𝐋, 𝐋, 𝐋]) ⊗ AxesDimensions(inv.([𝐓, 𝐋, 𝐋, 𝐋]))
```

or using [`nodims`](@ref) and scalar multiplication

```julia
𝐓 * nodims(2, 2)
# or dimensions(u"s") * nodims(2, 2)
```

for dimensionally homogeneous [`UnitfulTensor`](@ref)s.

## Defining functions of UnitfulTensors

Functions of [`UnitfulTensor`](@ref)s generally
act on the numerical values and dimensions independently
(with some exceptions, such as pivoted factorizations).
If you use a function from some package (or Julia Base)
that is not implemented in UnitfulTensors.jl, you can
implement it yourself.

For example, reversing a [`UnitfulTensor`](@ref)
along a specific dimension could be done like this:

```julia
using UnitfulTensors
import Base: reverse

function reverse(A::UnitfulTensor; kwargs...)
    dims = reverse(dimensions(A); kwargs...)
    vals = reverse(values(A); kwargs...)
    return UnitfulTensor(vals, dims)
end

function reverse(A::AxesDimensions; dims)
    d = dims
    scale, Adims... = dimsplat(A)
    newscale = scale * last(Adims[d]) # because AxisDimensions are normalized
    newdims = (Adims[1:d-1]..., reverse(Adims[d]), Adims[d+1:end]...)
    return AxesDimensions(newdims, newscale)
end

function reverse(A::AxisDimensions)
    dims = dimsvec(A)
    newdims = reverse(dims) ./ last(dims) # normalization
    return AxisDimensions(newdims)
end
```

Note that due to the tensor product structure of [`AxesDimensions`](@ref)
dealing with dimensions is usually much easier than with numerical values.

## Pitfalls

The types defined in this package store references to arrays of [`AbstractDimensions`](@ref).
Don't mutate them, or there will be surprises:

```jldoctest; setup = :(using UnitfulTensors)
julia> Z = UnitfulTensor([1u"S" 2u"S"
                          3u"S" 4u"S"])
       V = UnitfulTensor([5u"V", 6u"V"])
       I = Z * V
2-element UnitfulVector{Float64, SIDimensions, Vector{Float64}, AxesDimensions{1, SIDimensions}}:
 17.0 A
 39.0 A
 
julia> dimsvec(normdims(Z)[1])[2] = dimensions(u"V")
       I
2-element UnitfulVector{Float64, SIDimensions, Vector{Float64}, AxesDimensions{1, SIDimensions}}:
           17.0 A
 39.0 kg m^2 s^-3
```

Make a copy instead if necessary.

Currently dimensions are compared with `==` for performance reasons,
even though they are represented as floating-point numbers.
As a consequence, units with integer, half-integer, quarter-integer exponents work, while
units like s^(1/3) don't:

```julia-repl
julia> 1u"s^(2/3)" + 1u"s" / 1u"s^(1/3)"
ERROR: DimensionMismatch: dimensions (𝐓^11184811/16777216, 𝐓^5592405/8388608) do not match
```

This will be fixed when I figure out how to switch to `≈` without a significant performance penalty.

## Beyond SI

By default, UnitfulTensors.jl converts all quantities into SI units upon construction.
These units are provided by the `upreferred` function from [Unitful.jl](https://painterqubits.github.io/Unitful.jl).
It is possible to [override this behaviour](https://painterqubits.github.io/Unitful.jl/stable/conversion/#Fallback-promotion-rules)
and [define custom units](https://painterqubits.github.io/Unitful.jl/stable/newunits/)
by mechanisms described in the [Unitful.jl documentation](https://painterqubits.github.io/Unitful.jl).

Example:

```julia
import Unitful

const var"@uu_str" = Unitful.var"@u_str"
Unitful.@unit ħunit "ħ" ReducedPlanck Unitful.ħ false
Unitful.@unit eunit "e" ElementaryCharge Unitful.q false
Unitful.register(@__MODULE__)
_uT = ħunit/uu"meV"; _uL = uu"nm"; _uM = ħunit*_uT/_uL^2; _uI = eunit/_uT 
Unitful.preferunits(_uT, _uL, _uM, _uI)

using UnitfulTensors

σ = 1u"Ω^-1"; ω = 1u"THz2π"; E = 1u"eV"; vF = 1e6u"m/s"
σ, ω, E / vF

# output

(4108.235902227661 e^2 ħ^-1, 4.135667696923859 meV ħ^-1, 1.5192674478786259 ħ nm^-1)
```

However, this is not recommended unless really necessary. This feature has not been tested extensively
and it is better to stick to SI for consistency.

## Simplified UnitfulTensor construction

If you feel lazy to wrap all array literals in `UnitfulTensor(...)`, you can call `tensorize_literals()`.
After that, array literals containing [`UnitfulScalar`](@ref)s will return [`UnitfulTensor`](@ref)s
instead of `Array`s.

```jldoctest
using UnitfulTensors

tensorize_literals()

[1u"s" 2
 3     4u"s^-1"]
 
# output

2×2 UnitfulMatrix{Float64, SIDimensions, Matrix{Float64}, AxesDimensions{2, SIDimensions}}:
 1.0 s       2.0
   3.0  4.0 s^-1
```
 
!!! note
    This feature is experimental. It might break something and may be removed in the future.
