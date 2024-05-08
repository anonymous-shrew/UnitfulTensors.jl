for f in (:cond, :condskeel, :norm, :opnorm, :normalize, :normalize!)
    @eval @inline function $f(A::AbstractUnitfulScalarOrTensor; kwargs...)
        args = (A, )
        dims = $f(dimensions.(args)...; kwargs...)
        vals = $f(values.(args)...; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return T(vals, dims)
    end
    @eval @inline function $f(A::AbstractUnitfulScalarOrTensor, optarg::Real; kwargs...)
        args = (A, )
        dims = $f(dimensions.(args)..., optarg; kwargs...)
        vals = $f(values.(args)..., optarg; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return T(vals, dims)
    end
end

for f in (:condskeel, )
    @eval @inline function $f(A::AbstractUnitfulScalarOrTensor,
                              v::AbstractUnitfulScalarOrTensor; kwargs...)
        args = (A, v)
        dims = $f(dimensions.(args)...; kwargs...)
        vals = $f(values.(args)...; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return T(vals, dims)
    end
    @eval @inline function $f(A::AbstractUnitfulScalarOrTensor,
                              v::AbstractUnitfulScalarOrTensor, optarg::Real; kwargs...)
        args = (A, v)
        dims = $f(dimensions.(args)..., optarg; kwargs...)
        vals = $f(values.(args)..., optarg; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return T(vals, dims)
    end
end

rank(A::AbstractUnitfulScalarOrTensor; kwargs...) = rank(values(A); kwargs...)

for f in (:nullspace, )
    @eval @inline function $f(A::AbstractUnitfulScalarOrTensor; kwargs...)
        args = (A, )
        dims = $f(dimensions.(args)...; kwargs...)
        vals = $f(values.(args)...; kwargs...)
        n = size(vals, 2)
        T = promote_unitful($f, vals, dims ⊗ nodims(n), args...)
        return T(vals, dims ⊗ nodims(n))
    end
end

for f in (:norm, :opnorm)
    @eval function $f(A::AbstractDimensionsOrAxesDimensions, p::Real=2)
        if !ishomogeneous(A)
            throw(DimensionMismatch("$($f) requires a dimensionally homogeneous matrix"))
        end
        return dimscale(A)
    end
end

normalize(x::AbstractDimensionsOrAxesDimensions) = x / norm(x)
normalize(x::AbstractDimensionsOrAxesDimensions, p::Real) = x / norm(x, p)
normalize!(x::AbstractDimensionsOrAxesDimensions) = normalize(x)
normalize!(x::AbstractDimensionsOrAxesDimensions, p::Real) = normalize(x, p)

function cond(A::AbstractMatrixDimensions, p::Real=2)
    if !ishomogeneous(A)
        throw(DimensionMismatch("cond requires a homogeneous matrix"))
    else
        return one(eltype(A))
    end
end

function condskeel(A::AbstractMatrixDimensions, p::Real=Inf)
    if !ishomogeneous(A, 2)
        # TODO: fix the ordering of M and M^-1 in the Julia manual
        throw(DimensionMismatch("condskeel requires a matrix with dimensionally homogeneous rows"))
    else
        return one(eltype(A))
    end
end

function condskeel(A::AbstractMatrixDimensions, v::AbstractVectorDimensions, p::Real=Inf)
    if !ishomogeneous(A, 2)
        throw(DimensionMismatch("condskeel requires a matrix with dimensionally homogeneous rows"))
    elseif !ishomogeneous(v) || length(v) != size(A, 2)
        throw(DimensionMismatch("condskeel requires a dimensionally homogeneous vector \
                of length size(A, 2) == $(size(A, 2))"))
    else
        return one(eltype(A))
    end
end

function rank(A::AbstractDimensionsOrAxesDimensions; kwargs...)
    throw(ArgumentError("call rank directly on unitful matrices instead of their dimensions"))
end

function nullspace(A::AbstractMatrixDimensions; kwargs...)    
    return inv(normdims(A)[2])
end