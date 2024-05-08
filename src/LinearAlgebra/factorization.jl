# Could be <: AbstractUnitfulMatrix{T} if it's a factorization of a single matrix
# Pros: - a lot of methods for free
# Cons: - would require a separate AbstractUnitfulGeneralizedFactorization type for factorizations of two matrices
#       - inconsistent with the type hierarchy in stdlib
"""
    AbstractUnitfulFactorization{T}

Unitful version of `LinearAlgebra.Factorization{T}`.
"""
abstract type AbstractUnitfulFactorization{T} end

"""
    UnitfulFactorization{TV, TD, F<:Factorization{TV}, D<:AbstractMatrixDimensions{TD}} <:
                            AbstractUnitfulFactorization{UnitfulScalar{TV, TD}}

Factorization of an [`AbstractUnitfulMatrix`](@ref).
"""
struct UnitfulFactorization{TV, TD,
                            F<:Factorization{TV},
                            D<:AbstractMatrixDimensions{TD}} <:
                            AbstractUnitfulFactorization{UnitfulScalar{TV, TD}}
    numvals::F # values would conflict with eigen(A).values, etc.
    dims::D
end

"""
    UnitfulGeneralizedFactorization{TV, TD,
                            F<:Factorization{TV},
                            D<:Tuple{Vararg{AbstractMatrixDimensions{TD}}}
                            } <:
                            AbstractUnitfulFactorization{UnitfulScalar{TV, TD}}

Generalized factorization of two [`AbstractUnitfulMatrices`](@ref AbstractUnitfulMatrix).
"""
struct UnitfulGeneralizedFactorization{TV, TD,
                            F<:Factorization{TV},
                            D<:Tuple{Vararg{AbstractMatrixDimensions{TD}}}
                            } <:
                            AbstractUnitfulFactorization{UnitfulScalar{TV, TD}}
    numvals::F
    dims::D
end

values(A::AbstractUnitfulFactorization) = getfield(A, :numvals)
dimensions(A::AbstractUnitfulFactorization) = getfield(A, :dims)

# Can't use size(values(A)) because size(::Eigen) is undefined
size(A::AbstractUnitfulFactorization) = size(dimensions(A))

function iterate(F::AbstractUnitfulFactorization, args...)
    return iterate((getproperty(F, p) for p in _iterated_properties(F)), args...)
end

function getproperty(F::AbstractUnitfulFactorization, d::Symbol)
    if d ∈ _unitful_properties(F)
        Fvals = values(F); Fdims = dimensions(F)
        vals = getproperty(Fvals, d)
        dims = _factor_dimensions(Fvals, Fdims, d)
        T = promote_unitful(getproperty, vals, dims, F)
        return T(vals, dims)
    elseif d ∈ _unitless_properties(F)
        return getproperty(values(F), d)
    else
        return getfield(F, d)
    end
end

for f in (:_iterated_properties, :_displayed_properties, :_unitful_properties, :_unitless_properties)
    @eval $f(F::AbstractUnitfulFactorization) = $f(values(F))
end

_unitful_properties(F::Factorization) = _iterated_properties(F)
_displayed_properties(F::Factorization) = _iterated_properties(F)
_unitless_properties(::Factorization) = ()

_sym_uplo(F::Factorization) = sym_uplo(F.uplo)

issuccess(F::AbstractUnitfulFactorization; kwargs...) = issuccess(values(F); kwargs...)
