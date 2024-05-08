≈(xs::Vararg{AbstractUnitfulScalarOrTensor, 2}) = ≈(dimensions.(xs)...) && ≈(values.(xs)...)
≈(xs::Vararg{T, 2}) where T <: AbstractAxesDimensions = ≈(dimscale.(xs)...) && allmap(≈, normdims.(xs)...)
≈(xs::Vararg{T, 2}) where T <: AbstractAxisDimensions = allmap(≈, xs...)



function _mutability_warning(f)
    @warn "Dimensions of a UnitfulTensor are immutable. \
    Make sure you use the return value of $f instead of relying on argument mutation."
end

function _non_mutating(f)
    if f ∈ (:lmul!, :rmul!)
        return :*
    elseif f ∈ (:ldiv!, )
        return :\
    elseif f ∈ (:rdiv!, )
        return :/
    else
        throw(ArgumentError("Invalid argument $f in _non_mutating"))
    end
end

function _check_dims_pinv(A::AbstractVecOrMatDimensions)
    m, n = size(A, 1), size(A, 2)
    if m > n
        if !ishomogeneous(A, 1)
            throw(DimensionMismatch(
                    "matrix of an overdetermined system must have dimensionally homogeneous columns"))
        end
        if !ishomogeneous(A, 2)
            @warn "dimensionally inhomogeneous matrix of an overdetermined system must have full column rank. \
            $(@__MODULE__).jl doesn't check that."
        end
    elseif m < n
        if !ishomogeneous(A, 2)
            throw(DimensionMismatch(
                    "matrix of an underdetermined system must have dimensionally homogeneous rows"))
        end
        if !ishomogeneous(A, 1)
            @warn "dimensionally inhomogeneous matrix of an underdetermined system must have full row rank. \
            $(@__MODULE__).jl doesn't check that."
        end
    end
end

function _define_binary(f, T1, T2)
    return quote
        @inline function $f(A::$T1, B::$T2)
                    args = (A, B); kwargs = ()
                    dims = $f(dimensions.(args)...; kwargs...)
                    vals = $f(values.(args)...; kwargs...)
                    T = promote_unitful($f, vals, dims, args...)
                    return T(vals, dims)
                end
    end
end

function _define_binary_unitful_unitless(f, UT1, UT2)
    return quote
            $((_define_binary(f, T1, T2) for (T1, T2) in _unitful_unitless(UT1, UT2))...)
        end
end

_unitful_unitless(T1, T2) = ((T1, T2), (_unitless(T1), T2), (T1, _unitless(T2)))
_unitless(::Type{AbstractUnitfulScalar}) = Number
_unitless(::Type{AbstractUnitfulTensor}) = AbstractArray
_unitless(::Type{AbstractUnitfulVector}) = AbstractVector
_unitless(::Type{AbstractUnitfulMatrix}) = AbstractMatrix



for f in (:+, :-, :*, :\, :/, :dot, :cross, :kron, :inv, :pinv, 
          :lmul!, :ldiv!, :rmul!, :rdiv!,
          :sqrt)
    @eval @inline function $f(args::AbstractUnitfulScalarOrTensor...; kwargs...)
        dims = $f(dimensions.(args)...; kwargs...)
        vals = $f(values.(args)...; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return T(vals, dims)
    end
end

for f in (:+, :-, :*, :\, :/, :dot, :cross, :kron)
    for UT1 in (AbstractUnitfulScalar, AbstractUnitfulVector, AbstractUnitfulMatrix)
        for UT2 in (AbstractUnitfulScalar, AbstractUnitfulVector, AbstractUnitfulMatrix)
            @eval $(_define_binary_unitful_unitless(f, UT1, UT2))
        end
    end
end

for T1 in (AbstractUnitfulScalar, AbstractUnitfulMatrix)
    for T2 in (Integer, Rational, Real, Number)
        @eval function (^)(A::$T1, p::$T2)
            dims = dimensions(A)^p
            vals = values(A)^p
            T = promote_unitful(^, vals, dims, A)
            return T(vals, dims)
        end
    end
end

# Resolving ambiguities
for (T1, T2) in ((AdjointAbsVec{<:Number}, AbstractUnitfulVector),
                 (AdjOrTransAbsVec, AbstractUnitfulVector),
                 (AbstractUnitfulVector, AdjOrTransAbsVec),
                 (AdjointAbsVec, AbstractUnitfulMatrix),
                 (TransposeAbsVec, AbstractUnitfulMatrix))
    @eval $(_define_binary(:*, T1, T2))
end

for (T1, T2) in ((AdjointAbsVec, AbstractUnitfulMatrix), 
                 (TransposeAbsVec, AbstractUnitfulMatrix))
    @eval $(_define_binary(:/, T1, T2))
end

for f in (:+, :-)
    @eval $f(x::AbstractDimensionsOrAxesDimensions) = x
    @eval function $f(xs::Vararg{AbstractDimensionsOrAxesDimensions, 2})
        if xs[1] == xs[2]
            return xs[1]
        else
            throw(DimensionMismatch("dimensions $((xs[1], xs[2])) do not match"))
        end
    end
end

*(x::AbstractDimensions) = x
*(x::AbstractMatrixDimensions) = x

*(A::AbstractDimensions, B::AbstractAxisDimensions) = A ⊗ B
*(A::AbstractAxisDimensions, B::AbstractDimensions) = A ⊗ B

*(A::AbstractDimensions, B::AbstractAxesDimensions) = A ⊗ B
*(A::AbstractAxesDimensions, B::AbstractDimensions) = A ⊗ B

function *(A::AbstractMatrixDimensions, B::AbstractVecOrMatDimensions)
    scaleA, rowdimsA, coldimsA = dimsplat(A)
    scaleB, rowdimsB, maybe_coldimsB... = dimsplat(B)
    if !match(coldimsA, rowdimsB)
        throw(DimensionMismatch("column dimensions $coldimsA and row dimensions $rowdimsB are not multipliable"))
    end
    scale = scaleA * scaleB
    dims = (rowdimsA, maybe_coldimsB...)
    T = promote_dims(*, (dims, scale), A, B)
    return T(dims, scale)
end

function *(A::AbstractVectorDimensions, B::AbstractMatrixDimensions)
    scaleA, rowdimsA = dimsplat(A)
    scaleB, rowdimsB, coldimsB = dimsplat(B)
    if size(B, 1) != 1
        throw(DimensionMismatch("vector * matrix requires a single-row matrix"))
    end
    scale = scaleA * scaleB
    dims = (rowdimsA, coldimsB)
    T = promote_dims(*, (dims, scale), A, B)
    return T(dims, scale)
end

*(u::AdjointVectorDimensions, x::AbstractDimensions) = (x' * u')'
*(x::AbstractDimensions, u::AdjointVectorDimensions) = (u' * x')'
/(u::AdjointVectorDimensions, x::AbstractDimensions) = (x' \ u')'
\(x::AbstractDimensions, u::AdjointVectorDimensions) = (u' / x')'

*(u::AdjointVectorDimensions, v::AbstractVectorDimensions) = u' ⋅ v
*(u::AbstractVectorDimensions, v::AdjointVectorDimensions) = u ⊗ v'
*(u::AdjointVectorDimensions, v::AdjointVectorDimensions) = throw(MethodError(*, (u, v)))
*(u::AdjointVectorDimensions, A::AbstractMatrixDimensions) = (A' * u')'
/(u::AdjointVectorDimensions, A::AbstractMatrixDimensions) = (A' \ u')'
    
function \(A::AbstractMatrixDimensions, B::AbstractVecOrMatDimensions)
    scaleA, rowdimsA, coldimsA = dimsplat(A)
    scaleB, rowdimsB, maybe_coldimsB... = dimsplat(B)
    if rowdimsA != rowdimsB
        throw(DimensionMismatch("row dimensions of A and B do not match"))
    end
    _check_dims_pinv(A)
    scale = scaleA \ scaleB
    dims = (inv(coldimsA), maybe_coldimsB...)
    T = promote_dims(\, (dims, scale), A, B)
    return T(dims, scale)
end

function \(A::AbstractDimensions, B::AbstractAxesDimensions)
    scaleB, dims... = dimsplat(B)
    scale = A \ scaleB
    T = promote_dims(\, (dims, scale), A, B)
    return T(dims, scale)
end
    
function /(A::AbstractMatrixDimensions, B::AbstractMatrixDimensions)
    scaleA, rowdimsA, coldimsA = dimsplat(A)
    scaleB, rowdimsB, coldimsB = dimsplat(B)
    if coldimsA != coldimsB
        throw(DimensionMismatch("column dimensions of A and B do not match"))
    end
    _check_dims_pinv(B')
    scale = scaleA / scaleB
    dims = (rowdimsA, inv(rowdimsB))
    T = promote_dims(/, (dims, scale), A, B)
    return T(dims, scale)
end

/(A::AbstractAxesDimensions, B::AbstractDimensions) = B \ A

function inv(A::AbstractAxisDimensions)
    dims = inv.(A)
    T = promote_dims(inv, dims, A)
    return T(dims)
end

promote_dims(::typeof(inv), dims::Vector, A::AbstractAxisDimensions) = typeof(A)

function inv(A::AbstractMatrixDimensions)
    dims, scale = reverse(inv.(normdims(A))), inv(dimscale(A))
    T = promote_dims(inv, (dims, scale), A)
    return T(dims, scale)
end

pinv(A::AbstractMatrixDimensions; kwargs...) = (_check_dims_pinv(A); inv(A))

function pinv(A::AbstractVectorDimensions; kwargs...)
    _check_dims_pinv(A)
    dims, scale = inv.(normdims(A)), inv(dimscale(A))
    T = promote_dims(pinv, (dims, scale), A)
    return transpose(T(dims, scale))
end

pinv(A::AbstractDimensions; kwargs...) = inv(A)

for f in (:lmul!, :ldiv!, :rmul!, :rdiv!)
    g = _non_mutating(f)
    @eval function $f(A::AbstractDimensionsOrAxesDimensions, B::AbstractDimensionsOrAxesDimensions)
        _mutability_warning($f)      
        return $g(A, B)
    end
end

function ldiv!(Y::AbstractMatrixDimensions,
        A::Union{AbstractMatrixDimensions, AbstractDimensions}, B::AbstractMatrixDimensions)
    if Y != A \ B
        throw(DimensionMismatch("Dimensions of Y don't match A\\B"))
    end
    return Y
end

function dot(A::AbstractDimensionsOrAxesDimensions, B::AbstractDimensionsOrAxesDimensions)
    if ndims(A) != ndims(B)
        throw(DimensionMismatch("dot with unitful arguments of different ndims is not supported"))
    elseif !allmap(match, normdims(A), normdims(B))
        throw(DimensionMismatch("dot requires dimensionally reciprocal arguments (up to a scalar factor)"))
    end
    return dimscale(A) * dimscale(B)
end

function dot(x::AbstractDimensionsOrAxesDimensions,
             A::AbstractDimensionsOrAxesDimensions,
             y::AbstractDimensionsOrAxesDimensions)
    return dot(x, A*y)
end

function cross(A::AbstractVectorDimensions, B::AbstractVectorDimensions)
    if !(length(A) == length(B) == 3)
        throw(DimensionMismatch("cross product is only defined for vectors of length 3"))
    elseif !ishomogeneous(A) || !ishomogeneous(B)
        throw(DimensionMismatch("cross product of dimensionally inhomogeneous vectors is likely an error"))
    end
    return (dimscale(A) * dimscale(B)) ⊗ nodims(3)
end

function kron(A::AbstractAxisDimensions, B::AbstractAxisDimensions)
    dims = kron(dimsvec(A), dimsvec(B))
    T = promote_dims(kron, dims, A, B)
    return T(dims)
end

promote_dims(::typeof(kron), dims::Vector, args...) = AxisDimensions

_kron(A, B) = kron(A, B)
_kron(A, ::Nothing) = A
_kron(::Nothing, B) = B
_pad(t::Tuple, l::Integer) = (t..., ntuple(_ -> nothing, l - length(t))...)
_pad(ts::Tuple...) = _pad.(ts, maximum(length, ts))

function kron(A::AbstractVecOrMatDimensions, B::AbstractVecOrMatDimensions)
    scale, dims... = _kron.(_pad(dimsplat(A), dimsplat(B))...)
    T = promote_dims(kron, (dims, scale), A, B)
    return T(dims, scale)
end

function kron!(C::AbstractVecOrMatDimensions,
        A::AbstractDimensionsOrAxesDimensions, B::AbstractDimensionsOrAxesDimensions)
    if C != kron(A, B)
        throw(DimensionMismatch("kron!"))
    end
    return C
end

function kron!(C::AbstractUnitfulVecOrMat,
        A::AbstractUnitfulScalarOrTensor, B::AbstractUnitfulScalarOrTensor)
    kron!(dimensions.((C, A, B))...)
    kron!(values.((C, A, B))...)
    return C
end

kron(A::AdjointAxesDimensions, B::AdjointAxesDimensions) = kron(A', B')'

for T2 in (Integer, Real, Number)
    @eval function (^)(A::AbstractMatrixDimensions, p::$T2)
        scaleA, dims... = dimsplat(A)
        if !issquareable(A)
            throw(DimensionMismatch("only a squareable matrix can be raised to a power"))
        end
        scale = scaleA ^ p
        T = promote_dims(^, (dims, scale), A)
        return T(dims, scale)
    end
end

sqrt(A::AbstractMatrixDimensions) = A^(1//2)