# How indexing into AbstractAxesDimensions works:
#
# getindex(A, I1, I2, ...) invokes all the to_indices/checkbounds machinery of Julia Base
# until Base._unsafe_getindex is reached.
#
# Base._unsafe_getindex is overloaded in this file to take into account
# the tensor product structure of AbstractAxesDimensions.
#
# Each index acts on a certain dimension or dimensions of A. (dimensions in the mathematical sense)
# Their number is the "output dimension" of an index and is obtained with Base.index_ndims.
# E. g., "output dimension" = 1 for numbers and arrays thereof;
#        "output dimension" = N for CartesianIndex{N} and CartesianIndices{N}.
#
# Tensor product structure allows us to factorize A
# and process each index separately with _get_single_index:
# A[I1, I2, ...] = dimscale(A) * A1[I1] ⊗ A2[I2] ⊗ ..
#
# Since A1[I1], A2[I2], ... also must have a tensor product structure,
# we compute their first row/column/etc. with _getAxesDims/_getAxisDims
# and check that the remaining elements obey the tensor product structure.
# This check is done with _setindex_rest!, which calls Base._unsafe_getindex!,
# which calls setindex! on every element of A1[I1], etc., which throws an error
# whenever the tensor product structure is violated.

"""
`_N_but_one(4) == ((2, 3, 4), (1, 3, 4), (1, 2, 4), (1, 2, 3))`
"""
_N_but_one(N, k) = (1:(k-1)..., (k+1):N...)
_N_but_one(N) = _N_but_one.(N, (1:N...,))

# LogicalIndex is a non-indexable AbstractVector, which confuses Iterators.pairs (and me)
_pairs(I::AbstractVector) = zip(eachindex(I), I)



IndexStyle(A::AbstractAxisDimensions) = IndexLinear()
IndexStyle(A::AbstractAxesDimensions) = IndexCartesian()
IndexStyle(A::AbstractUnitfulTensor) = IndexCartesian()

size(A::AbstractAxisDimensions) = (length(dimsvec(A)),) 
size(A::AbstractAxesDimensions) = length.(normdims(A))
size(A::AbstractUnitfulTensor) = size(values(A))



function getindex(A::AbstractAxisDimensions, i::Int)
    return dimsvec(A)[i]
end

function getindex(A::AbstractAxisDimensions, I::AbstractVector)
    scale = isempty(I) ? one(eltype(A)) : A[first(I)]
    dims = similar(A, axes(I))
    for (i, Ii) in _pairs(I)
        dims[i] = A[Ii] / scale
    end
    T = promote_dims(getindex, ((dims, ), scale), A)
    return T((dims, ), scale; normalize = false)
end



function getindex(A::AbstractAxesDimensions{N}, I::Vararg{Int, N}) where N
    scale, dims... = dimsplat(A)
    return scale * prod(getindex.(dims, I); init = one(eltype(A)))
end

function setindex!(A::AbstractAxesDimensions, x, I...)
    x == A[I...] || throw(ArgumentError("Changing the dimensions of AxesDimensions with setindex! is not supported"))
    return x
end

@inline function _unsafe_getindex(::IndexStyle, A::AbstractAxesDimensions, I::Union{Real, AbstractArray}...)
    inds_out_dims = length.(index_ndims.(I))
    factors, scale = tensor_factorize(A, inds_out_dims)
    # map is faster than broadcasting
    return scale * tensor_product(map(_get_single_index, factors, I)...; init = one(eltype(A)))
end

# _get_single_index and related functions can be used to convert
# AbstractArray{<:AbstractDimensions} into AbstractAxesDimensions,
# so they should accept AbstractDimensionsArray, not only AbstractAxesDimensions
const AbstractDimensionsArray = AbstractArray{<:AbstractDimensions}

@inline function _get_single_index(A::AbstractDimensionsArray, I::Union{Real, AbstractArray})
    isempty(I) && return nodims(size(I)...)
    scale = A[first(I)]
    dims = _getAxesDims(A, I)
    T = promote_dims(_get_single_index, (dims, scale), A)
    dest = T(dims, scale)
    return _setindex_rest!(dest, A, I)
end

_setindex_rest!(dest::AbstractAxesDimensions{0}, A::AbstractDimensionsArray,
    I::Union{Real, AbstractArray{<:Any, 0}}) = dest
_setindex_rest!(dest::AbstractVectorDimensions, A::AbstractDimensionsArray, I::AbstractVector) = dest
_setindex_rest!(dest::AbstractAxesDimensions, A::AbstractDimensionsArray,
    I::Union{Real, AbstractArray}) = _unsafe_getindex!(dest, A, I)

@inline function _getAxesDims(A::AbstractDimensionsArray, I::Union{Real, AbstractArray})
    @inline f(dims) = _getAxisDims(A, first(eachslice(I, dims = dims)))
    return map(f, _N_but_one(ndims(I)))
end

function _getAxisDims(A::AbstractDimensionsArray, I::AbstractVector)
    scale = isempty(I) ? one(eltype(A)) : _getindex_normdims(A, first(I))
    shape = index_shape(I)
    newdims = similar(A, shape)
    for (i, Ii) in _pairs(I)
         newdims[i] = _getindex_normdims(A, Ii) / scale
    end
    T = promote_dims(_getAxisDims, newdims, A)
    return T(newdims; normalize = false)
end

promote_dims(::typeof(_getAxisDims), dims, args...) = AxisDimensions

@inline _prod(itr; init) = isempty(itr) ? init : prod(itr) # this cuts 4-6 ns in some benchmarks

function _getindex_normdims(A::AbstractAxesDimensions, I)
    dims = normdims(A)
    inds = to_indices(A, (I,))
    return _prod(getindex.(dims, inds); init = one(eltype(A)))
end

_getindex_normdims(A::AbstractDimensionsArray, I) = A[I]



#### Special cases ####

## Omitted singleton dimensions and trailing 1's ##

const _singleton_dim = AxisDimensions([NoDims])

function reshape(A::AbstractAxesDimensions{M}, ndims::Val{N}) where {M, N}
    scale, dims... = dimsplat(A)
    if M > N
        all(length.(dims[(N+1):end]) .== 1) || throw(ArgumentError(
                "Only trailing dimensions of length 1 can be removed by reshape(::AxesDimensions, ::Val)"))
        newdims = dims[1:N]
    else
        newdims = (dims..., ntuple(_ -> _singleton_dim, N - M)...)
    end
    T = promote_dims(reshape, (newdims, scale), A)
    return T(newdims, scale)
end

# Overloading Base._maybe_reshape(::IndexCartesian, A::AbstractVector, I...) = A
_maybe_reshape(::IndexCartesian, A::AbstractVectorDimensions, I...) = __maybe_reshape(A, index_ndims(I...))

## Logical indexing ##

first(I::LogicalIndex) = first(i for i in I)
_getAxesDims(A::AbstractAxesDimensions, I::LogicalIndex) = (_getAxisDims(A, I), )

## Optimization for Colon ##

function getindex(A::AbstractAxisDimensions, I::Slice)
    dims = (copy(dimsvec(A)),)
    T = promote_dims(getindex, dims, A)
    return T(dims; normalize = false)
end

_get_single_index(A::AbstractVectorDimensions, I::Slice) = normdims(A)[1][I]

################################################################################

@inline function getindex(A::UnitfulTensor, inds...)
    dims = getindex(dimensions(A), inds...)
    vals = getindex(values(A), inds...)
    T = promote_unitful(getindex, vals, dims, A)
    return T(vals, dims)
end

@inline function setindex!(A::UnitfulTensor, v, inds...)
    dims = setindex!(dimensions(A), dimensions(v), inds...)
    vals = setindex!(values(A), values(v), inds...)
    return v
end