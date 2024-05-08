for f in (:adjoint, :adjoint!, :transpose, :transpose!)
    @eval @inline function $f(args::AbstractUnitfulScalarOrTensor...; kwargs...)
        dims = $f(dimensions.(args)...; kwargs...)
        vals = $f(values.(args)...; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return T(vals, dims)
    end
end

"""
    AdjointAxesDimensions{T<:AbstractDimensions, P<:AbstractVecOrMatDimensions{T}} <:
                          AbstractMatrixDimensions{T}

Like `LinearAlgebra.Adjoint` and `LinearAlgebra.Transpose`, but for [`AbstractAxesDimensions`](@ref).
"""
struct AdjointAxesDimensions{T<:AbstractDimensions, P<:AbstractVecOrMatDimensions{T}} <:
                             AbstractMatrixDimensions{T}
    parent::P
end

"""
    AdjointVectorDimensions{T<:AbstractDimensions}

Alias for [`AdjointAxesDimensions{T, P} where P <: AbstractVectorDimensions{T}`](@ref).
"""
const AdjointVectorDimensions{T<:AbstractDimensions} = AdjointAxesDimensions{T, P} where
                                                       P <: AbstractVectorDimensions{T}
                                                       
"""
    AdjointMatrixDimensions{T<:AbstractDimensions}

Alias for [`AdjointAxesDimensions{T, P} where P <: AbstractMatrixDimensions{T}`](@ref).
"""
const AdjointMatrixDimensions{T<:AbstractDimensions} = AdjointAxesDimensions{T, P} where
                                                       P <: AbstractMatrixDimensions{T}


parent(x::AdjointAxesDimensions) = x.parent

normdims(x::AdjointMatrixDimensions) = reverse(normdims(parent(x)))
normdims(x::AdjointVectorDimensions) = (_singleton_dim, normdims(parent(x))...)
dimscale(x::AdjointAxesDimensions) = dimscale(parent(x))

dimensions(A::Union{Adjoint, Transpose}) = nodims(size(parent(A))...)'



adjoint(A::AbstractVecOrMatDimensions) = AdjointAxesDimensions(A)
adjoint(x::AbstractDimensions) = x
adjoint(A::AdjointAxesDimensions) = parent(A)

function adjoint!(dest::AbstractAxesDimensions, src::AbstractAxesDimensions)
    if dest != src
        throw(DimensionMismatch("dimensions $dest and $src don't match"))
    else
        return dest
    end
end

transpose(A::AbstractDimensionsOrAxesDimensions) = adjoint(A)
transpose!(A::AbstractDimensionsOrAxesDimensions, B::AbstractDimensionsOrAxesDimensions) = adjoint!(A, B)
