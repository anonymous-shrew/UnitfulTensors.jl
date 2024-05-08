const AbstractAxesDimensionsLike = Union{AbstractAxesDimensions, AbstractAxisDimensions, AbstractDimensions}

"""
    tensor_product(xs::AbstractAxesDimensionsLike...; init = NoDims)
Compute the tensor product of [`AbstractAxesDimensions`](@ref),
[`AbstractAxisDimensions`](@ref), or [`AbstractDimensions`](@ref).

If the argument list may be empty, the init keyword must be supplied
([`NoDims`](@ref) if one works with [`SIDimensions`](@ref)).
`⊗` can be used as a shorthand.

# Examples:
```jldoctest
julia> 𝐋 ⊗ AxisDimensions([NoDims, 𝐌]) ⊗ AxisDimensions([NoDims, 𝐓])
2×2 AxesDimensions{2, SIDimensions}:
   𝐋    𝐋 𝐓
 𝐋 𝐌  𝐋 𝐌 𝐓
```
"""
function tensor_product(xs::AbstractAxesDimensionsLike...; init = NoDims)
    dims = flatten(map(normdims, xs))
    scale = prod(dimscale.(xs); init = init)
    T = promote_dims(tensor_product, (dims, scale), xs...)
    return T(dims, scale)
end

tensor_product() = throw(ArgumentError("tensor_product of zero arguments is undefined without the init keyword"))

const ⊗ = tensor_product

"""
    tensor_factorize(A::AbstractAxesDimensions, factor_ndims::Tuple{Vararg{Int}})
Represent `A` as `dimscale(A) ⊗ factors[1] ⊗ factors[2] ⊗ ...`,
where `dimscale(factors[i]) == NoDims` and `ndims.(factors) == factor_ndims`.

Returns `(factors, dimscale(A))`.
"""
@inline function tensor_factorize(A::AbstractAxesDimensions, factor_ndims::Tuple{Vararg{Int}})
    dimgroups = partition(normdims(A), factor_ndims)
    T = promote_dims(tensor_factorize, dimgroups, A)
    factors = T.(dimgroups)
    scale = dimscale(A)
    return factors, scale
end

"""
Partition a tuple into subtuples of given lengths.
"""
partition(t::Tuple, lengths::Tuple{}) = ()

@inline function partition(t::Tuple, lengths::Tuple{Vararg{Integer}})
    head, rest = split(t, Val(lengths[1]))
    return (head, partition(rest, tail(lengths))...)
end
