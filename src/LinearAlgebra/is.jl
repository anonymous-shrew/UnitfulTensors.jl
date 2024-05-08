for f in (:issymmetric, :ishermitian,
          :issquareable, :isendomorphic)
    @eval $f(A::AbstractUnitfulScalarOrTensor) = $f(dimensions(A)) && $f(values(A))
end

for f in (:isposdef, :isposdef!)
    @eval $f(A::AbstractUnitfulScalarOrTensor) = issymmetric(dimensions(A)) && $f(values(A))
end

for f in (:ishomogeneous, )
    @eval $f(A::AbstractUnitfulScalarOrTensor) = $f(dimensions(A)) && $f(values(A))
    @eval $f(A::AbstractUnitfulScalarOrTensor, d::Integer) = $f(dimensions(A), d) && $f(values(A), d)
end

for f in (:isdiag, )
    @eval $f(A::AbstractUnitfulScalarOrTensor) = $f(values(A))
end

for f in (:istril, :istriu)
    @eval $f(A::AbstractUnitfulScalarOrTensor) = $f(values(A))
    @eval $f(A::AbstractUnitfulScalarOrTensor, k::Integer) = $f(values(A), k)
end

issymmetric(A::AbstractMatrixDimensions) = normdims(A)[1] == normdims(A)[2]
issymmetric(A::AbstractDimensions) = true
ishermitian(A::AbstractDimensionsOrAxesDimensions) = issymmetric(A)

"""
    issquareable(A::AbstractUnitfulMatrix) -> Bool
    issquareable(A::AbstractMatrixDimensions) -> Bool
    issquareable(A::AbstractMatrix{<:Number}) -> Bool

Test whether a matrix can be squared.

Certain functions, such as `eigen` and `tr`, are defined only for squareable matrices.
"""
issquareable(A::AbstractMatrixDimensions) = match(normdims(A)[1], normdims(A)[2])

"""
    isendomorphic(A::AbstractUnitfulMatrix) -> Bool
    isendomorphic(A::AbstractMatrixDimensions) -> Bool
    isendomorphic(A::AbstractMatrix{<:Number}) -> Bool

Test whether a matrix is endomorphic (can represent a linear map from a vector space to itself).

Transcendental functions, such as `exp` and `log`, are defined only for endomorphic matrices.
"""
isendomorphic(A::AbstractMatrixDimensions) = (dimscale(A) == one(eltype(A))) && issquareable(A)

"""
    ishomogeneous(A::AbstractUnitfulTensor[, d::Integer]) -> Bool
    ishomogeneous(A::AbstractAxesDimensions[, d::Integer]) -> Bool
    ishomogeneous(A::AbstractArray{<:Number}[, d::Integer]) -> Bool
    ishomogeneous(A::AbstractAxisDimensions) -> Bool

Test whether an array is dimensionally homogeneous (all of its elements have the same physical dimensions).

`d` can be provided to test homogeneity along a specific (mathematical) dimension.

Most of linear algebra can be applied to dimensionally homogeneous matrices.
Notable functions that require homogeneity include `svd` and `norm`.
Least-squares solution of an overdetermined system with `\\` or `pinv`
requires homogeneity along the first dimension; otherwise the sum of squares is undefined.
"""
ishomogeneous(A::AbstractAxesDimensions) = all(ishomogeneous(axis) for axis in normdims(A))
ishomogeneous(A::AbstractAxesDimensions, d::Integer) = d > ndims(A) || ishomogeneous(normdims(A)[d])
ishomogeneous(A::AbstractAxisDimensions) = all(dims == one(eltype(A)) for dims in A)

for f in (:issquareable, :isendomorphic)
    @eval $f(A::Union{AbstractMatrix{<:Number}, Number}) = size(A, 1) == size(A, 2)
end

for f in (:ishomogeneous, )
    @eval $f(A::Union{AbstractArray{<:Number}, Number}) = true
    @eval $f(A::Union{AbstractArray{<:Number}, Number}, d::Integer) = true
end