function _pivots(f)
    if f ∈ (:bunchkaufman, :bunchkaufman!)
        return (:Bool, ) # bunchkaufman(A, rook::Bool=false; check = true)
    elseif f ∈ (:cholesky, :cholesky!)
        return (:NoPivot, :RowMaximum)
    elseif f ∈ (:lu, :lu!)
        return (:NoPivot, :RowMaximum, :RowNonZero)
    elseif f ∈ (:qr, :qr!)
        return (:NoPivot, :ColumnNorm)
    else
        return ()
    end
end

for f in (:lu, :lu!, 
          :bunchkaufman, :bunchkaufman!, :cholesky, :cholesky!, :ldlt, :ldlt!,
          :eigen, :eigen!, :hessenberg, :hessenberg!, :schur, :schur!, :svd, :svd!,
          :lq, :lq!, :qr, :qr!)
    @eval @inline function $f(A::AbstractUnitfulMatrix; kwargs...)
        vals = values(A); dims = dimensions(A)
        _check_dims($f, dims)
        return UnitfulFactorization($f(vals; kwargs...), dims)
    end
    # Required to resolve method ambiguities
    for pivot in _pivots(f)
        @eval @inline function $f(A::AbstractUnitfulMatrix, pivot::$pivot; kwargs...)
            vals = values(A); dims = dimensions(A)
            _check_dims($f, dims)
            return UnitfulFactorization($f(vals, pivot; kwargs...), dims)
        end
    end
end

# Generalized factorizations
for f in (:eigen, :eigen!, :schur, :schur!, :svd, :svd!)
    @eval @inline function $f(A::AbstractUnitfulMatrix, B::AbstractUnitfulMatrix; kwargs...)
        vals = values.((A, B)); dims = dimensions.((A, B))
        _check_dims($f, dims...)
        return UnitfulGeneralizedFactorization($f(vals...; kwargs...), dims)
    end
end

for f in (:eigvals, :eigvals!, :svdvals, :svdvals!, :eigvecs, :eigmax, :eigmin)
    @eval @inline function $f(args::AbstractUnitfulMatrix...; kwargs...)
        dims = $f(dimensions.(args)...; kwargs...)
        vals = $f(values.(args)...; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return T(vals, dims)
    end
end

for f in (:svdvals, :svdvals!)
    @eval @inline function $f(args::Vararg{AbstractUnitfulMatrix, 2}; kwargs...)
        _check_dims($f, dimensions.(args)...)
        return $f(values.(args)...; kwargs...)
    end
end



for f in (:eigvals, :eigvals!, :svdvals, :svdvals!)
    @eval function $f(A::AbstractMatrixDimensions; kwargs...)
        _check_dims($f, A)
        N = min(size(A)...)
        scale = dimscale(A)
        return scale ⊗ nodims(N)
    end
end

for f in (:eigvecs, )
    @eval function $f(A::AbstractMatrixDimensions; kwargs...)
        _check_dims($f, A)
        N = size(A, 1)
        rowdims = normdims(A)[1]
        return rowdims ⊗ nodims(N)
    end
end

for f in (:eigmax, :eigmin)
    @eval function $f(A::AbstractMatrixDimensions; kwargs...)
        _check_dims($f, A)
        return dimscale(A)
    end
end

for f in (:eigvals, :eigvals!)
    @eval function $f(A::AbstractMatrixDimensions, B::AbstractMatrixDimensions; kwargs...)
        _check_dims($f, A, B)
        N = size(A, 1)
        scaleA, scaleB = dimscale.((A, B))
        scale = scaleA/scaleB
        return scale ⊗ nodims(N)
    end
end

for f in (:eigvecs, )
    @eval function $f(A::AbstractMatrixDimensions, B::AbstractMatrixDimensions; kwargs...)
        _check_dims($f, A, B)
        N = size(A, 1)
        coldims = normdims(A)[2]
        return inv(coldims) ⊗ nodims(N)
    end
end



# By default, _unitful_properties(F) == _displayed_properties(F) == _iterated_properties(F)
# and _unitless_properties(F) == ()
_iterated_properties(F::BunchKaufman) = (:D, _sym_uplo(F), :p)
 _unitful_properties(F::BunchKaufman) = (:D, _sym_uplo(F), :P)
_unitless_properties(F::BunchKaufman) = (:p, :uplo, )
 _iterated_properties(F::Cholesky) = (:L, :U)
_displayed_properties(F::Cholesky) = (_sym_uplo(F), )
 _unitless_properties(F::Cholesky) = (:uplo,)
 _iterated_properties(F::CholeskyPivoted) = (:L, :U)
_displayed_properties(F::CholeskyPivoted) = (_sym_uplo(F), :p)
  _unitful_properties(F::CholeskyPivoted) = (:L, :U, :P)
 _unitless_properties(F::CholeskyPivoted) = (:p, :uplo)
_iterated_properties(F::Union{Eigen, GeneralizedEigen})  = (:values, :vectors)
 _iterated_properties(F::Hessenberg) = (:Q, :H, :μ)
_displayed_properties(F::Hessenberg) = (:μ, :Q, :H)
  _unitful_properties(F::Hessenberg) = (F.μ isa Bool) ? (:Q, :H) : (:Q, :H, :μ)
 _unitless_properties(F::Hessenberg) = (F.μ isa Bool) ? (:μ, ) : ()
 _iterated_properties(F::LDLt) = ()
_displayed_properties(F::LDLt) = (:L, :D)
  _unitful_properties(F::LDLt) = (:L, :D, :Lt, :d)
_iterated_properties(F::LQ) = (:L, :Q)
 _iterated_properties(F::LU) = (:L, :U, :p)
_displayed_properties(F::LU) = (:L, :U)
 _unitful_properties(F::LU) = (:L, :U, :P)
_unitless_properties(F::LU) = (:p, )
_iterated_properties(F::Union{QR, QRCompactWY}) = (:Q, :R)
_iterated_properties(F::QRPivoted) = (:Q, :R, :p)
 _unitful_properties(F::QRPivoted) = (:Q, :R, :P)
_unitless_properties(F::QRPivoted) = (:p, )
_iterated_properties(F::Schur) = (:T, :Z, :values)
 _unitful_properties(F::Schur) = (:Schur, :T, :vectors, :Z, :values)
_iterated_properties(F::GeneralizedSchur) = (:S, :T, :Q, :Z, :α, :β)
_unitful_properties(F::GeneralizedSchur)  = (:S, :T, :left, :Q, :right, :Z, :values, :α, :β)
 _iterated_properties(F::SVD) = (:U, :S, :V)
_displayed_properties(F::SVD) = (:U, :S, :Vt)
  _unitful_properties(F::SVD) = (:U, :S, :V, :Vt)
_iterated_properties(F::GeneralizedSVD) = (:U, :V, :Q, :D1, :D2, :R0)

_property_descriptions(F::CholeskyPivoted) = (; _sym_uplo(F) => "$(F.uplo) factor with rank $(rank(F))")
_property_descriptions(F::Union{Eigen, GeneralizedEigen}) = (values = "values", vectors = "vectors")
_property_descriptions(F::Schur) = (values = "eigenvalues", )
_property_descriptions(F::GeneralizedSchur) = (α = "α", β = "β")
_property_descriptions(F::SVD) = (S = "singular values", )

function _show_property(io::IO, mime::MIME{Symbol("text/plain")},
        F::UnitfulFactorization{<:Any,<:Any,<:Hessenberg}, ::Val{:μ})    
    if !iszero(F.μ)
        print("\nwith shift μI for μ = ", F.μ)
    end
end

_issuccess(F::BunchKaufman) = issuccess(F)
_issuccess(F::Union{Cholesky, CholeskyPivoted}) = issuccess(F)
_issuccess(F::LU) = F.info >= 0



function _left_permutation_matrix(A::AbstractAxesDimensions, p::AbstractVector)
    dims = normdims(A)[1]
    return dims[p] ⊗ inv(dims)
end

function _right_permutation_matrix(A::AbstractAxesDimensions, p::AbstractVector)
    dims = normdims(A)[2]
    return inv(dims) ⊗ dims[p]
end

function _throw_factor_dimensions(F, d)
    throw(ArgumentError("Error in _factor_dimensions(...):\
            type $(nameof(typeof(F))) doesn't have field $d or it is unitless."))
end

#### LU ####

function _check_dims(::Union{typeof.((lu, lu!))...}, dims::AbstractMatrixDimensions)
    return nothing
end

function _factor_dimensions(F::LU, A::AbstractAxesDimensions, d::Symbol)
    p = F.p
    K = min(size(A)...)
    scale, rowdims, coldims = dimsplat(A)
    if d === :L
        return sqrt(scale) ⊗ rowdims[p] ⊗ nodims(K)
    elseif d === :U
        return sqrt(scale) ⊗ nodims(K) ⊗ coldims
    elseif d === :P
        return _left_permutation_matrix(A, p)
    else
        _throw_factor_dimensions(F, d)
    end
end

#### Cholesky and friends ####

function _check_dims(f::Union{typeof.((bunchkaufman, bunchkaufman!, cholesky, cholesky!, ldlt, ldlt!))...},
        dims::AbstractMatrixDimensions)
    if !issymmetric(dims)
        throw(DimensionMismatch("$f requires a dimensionally symmetric matrix"))
    end
end

function _factor_dimensions(F::Cholesky, A::AbstractAxesDimensions, d::Symbol)
    N = size(A, 1)
    scale, rowdims, coldims = dimsplat(A)
    if d === :L
        return sqrt(scale) ⊗ rowdims ⊗ nodims(N)
    elseif d === :U
        return sqrt(scale) ⊗ nodims(N) ⊗ coldims
    else
        _throw_factor_dimensions(F, d)
    end
end

function _factor_dimensions(F::CholeskyPivoted, A::AbstractAxesDimensions, d::Symbol)
    p = F.p
    N = size(A, 1)
    scale, rowdims, coldims = dimsplat(A)
    if d === :L
        return sqrt(scale) ⊗ rowdims[p] ⊗ nodims(N)
    elseif d === :U
        return sqrt(scale) ⊗ nodims(N) ⊗ coldims[p]
    elseif d === :P
        return _left_permutation_matrix(A, p)
    else
        _throw_factor_dimensions(F, d)
    end
end

function _factor_dimensions(F::LDLt, A::AbstractAxesDimensions, d::Symbol)
    scale, rowdims, coldims = dimsplat(A)
    if d === :L
        return rowdims ⊗ inv(rowdims)
    elseif d === :Lt
        return inv(coldims) ⊗ coldims
    elseif d === :D
        return A
    elseif d === :d
        return diag(A)
    else
        _throw_factor_dimensions(F, d)
    end
end

function _factor_dimensions(F::BunchKaufman, A::AbstractAxesDimensions, d::Symbol)
    p = F.p
    scale, rowdims, coldims = dimsplat(A)
    if d === _sym_uplo(F)
        return rowdims[p] ⊗ inv(rowdims)[p]
    elseif d === :D
        return A[p, p]
    elseif d === :P
        return _left_permutation_matrix(A, p)
    else
        _throw_factor_dimensions(F, d)
    end
end

#### Eigen and friends ####

function _check_dims(f::Union{typeof.((eigen, eigen!, hessenberg, hessenberg!, schur, schur!,
                                       eigmax, eigmin, eigvals, eigvals!, eigvecs))...},
                     dims::AbstractMatrixDimensions)
    if !issquareable(dims)
        throw(DimensionMismatch("$f requires a squareable matrix"))
    end
end

function _check_dims(f::Union{typeof.((eigen, eigen!, schur, schur!,
                                       eigvals, eigvals!, eigvecs))...},
                     dimsA::AbstractMatrixDimensions, dimsB::AbstractMatrixDimensions)
    checksquare(dimsA); checksquare(dimsB)
    if normdims(dimsA) != normdims(dimsB)
        throw(DimensionMismatch("$f requires matrices of the same physical dimensions \
                (up to an overall scalar factor)"))
    end
end

function _factor_dimensions(F::Eigen, A::AbstractAxesDimensions, d::Symbol)
    N = size(A, 1)
    scale, rowdims, coldims = dimsplat(A)
    if d === :values
        return scale ⊗ nodims(N)
    elseif d === :vectors
        return rowdims ⊗ nodims(N)
    else
        _throw_factor_dimensions(F, d)
    end
end

function _factor_dimensions(F::GeneralizedEigen, AB::NTuple{2, AbstractAxesDimensions}, d::Symbol)
    A, B = AB
    N = size(A, 1)
    scaleA, rowdims, coldims = dimsplat(A)
    scaleB = dimscale(B)
    if d === :values
        scale = scaleA/scaleB
        return scale ⊗ nodims(N)
    elseif d === :vectors
        return inv(coldims) ⊗ nodims(N)
    else
        _throw_factor_dimensions(F, d)
    end
end

function _factor_dimensions(F::Hessenberg, A::AbstractAxesDimensions, d::Symbol)
    N = size(A, 1)
    scale, rowdims, coldims = dimsplat(A)
    if d === :μ
        return scale ⊗ nodims(N)
    elseif d === :Q
        return rowdims ⊗ nodims(N)
    elseif d === :H
        return scale ⊗ nodims(N, N)
    else
        _throw_factor_dimensions(F, d)
    end
end

function _factor_dimensions(F::Schur, A::AbstractAxesDimensions, d::Symbol)
    N = size(A, 1)
    scale, rowdims, coldims = dimsplat(A)
    if d === :values
        return scale ⊗ nodims(N)
    elseif d ∈ (:vectors, :Z)
        return rowdims ⊗ nodims(N)
    elseif d ∈ (:Schur, :T)
        return scale ⊗ nodims(N, N)
    else
        _throw_factor_dimensions(F, d)
    end
end

function _factor_dimensions(F::GeneralizedSchur, AB::NTuple{2, AbstractAxesDimensions}, d::Symbol)
    A, B = AB
    N = size(A, 1)
    scaleA, rowdims, coldims = dimsplat(A)
    scaleB = dimscale(B)
    if d === :values
        scale = scaleA/scaleB
        return scale ⊗ nodims(N)
    elseif d === :α
        return scaleA ⊗ nodims(N)
    elseif d === :β
        return scaleB ⊗ nodims(N)
    elseif d ∈ (:left, :Q)
        return rowdims ⊗ nodims(N)
    elseif d ∈ (:right, :Z)
        return coldims ⊗ nodims(N)
    elseif d === :S
        return scaleA ⊗ nodims(N, N)
    elseif d === :T
        return scaleB ⊗ nodims(N, N)
    else
        _throw_factor_dimensions(F, d)
    end
end

function _check_dims(f::Union{typeof.((svd, svd!, svdvals, svdvals!))...},
                     dims::AbstractMatrixDimensions)
    if !ishomogeneous(dims)
        throw(DimensionMismatch("$f requires a dimensionally homogeneous matrix"))
    end
end

function _check_dims(f::Union{typeof.((svd, svd!, svdvals, svdvals!))...},
                     dimsA::AbstractMatrixDimensions, dimsB::AbstractMatrixDimensions)
    if size(dimsA, 2) != size(dimsB, 2) || dimscale(dimsA) != dimscale(dimsB) || 
        !ishomogeneous(dimsA) || !ishomogeneous(dimsB)
        throw(DimensionMismatch("$f requires dimensionally homogeneous matrices \
                with the same number of columns and physical dimensions"))
    end
end
    
function _factor_dimensions(F::SVD, A::AbstractAxesDimensions, d::Symbol)
    MU, NU = size(F.U)
    MVt, NVt = size(F.Vt)
    NS = length(F.S)
    scale = dimscale(A)
    if d === :U
        return nodims(MU, NU)
    elseif d === :Vt
        return nodims(MVt, NVt)
    elseif d === :V
        return nodims(NVt, MVt)
    elseif d === :S
        return scale ⊗ nodims(NS)
    else
        _throw_factor_dimensions(F, d)
    end
end

function _factor_dimensions(F::GeneralizedSVD, AB::NTuple{2, AbstractAxesDimensions}, d::Symbol)
    A, B = AB
    M, N = size(A)
    P = size(B, 1)
    K, L = F.k, F.l
    scale = dimscale(A)
    if d === :U
        return nodims(M, M)
    elseif d === :V
        return nodims(P, P)
    elseif d === :Q
        return nodims(N, N)
    elseif d === :D1
        return nodims(M, K + L)
    elseif d === :D2
        return nodims(P, K + L)
    elseif d === :R0
        return scale ⊗ nodims(K + L, N)
    else
        _throw_factor_dimensions(F, d)
    end
end

#### QR and friends ####

function _check_dims(f::Union{typeof.((lq, lq!))...}, dims::AbstractMatrixDimensions)
    if !ishomogeneous(dims, 2)
        throw(DimensionMismatch("$f requires a matrix with dimensionally homogeneous rows"))
    end
end

function _check_dims(f::Union{typeof.((qr, qr!))...}, dims::AbstractMatrixDimensions)
    if !ishomogeneous(dims, 1)
        throw(DimensionMismatch("$f requires a matrix with dimensionally homogeneous columns"))
    end
end

function _factor_dimensions(F::LQ, A::AbstractAxesDimensions, d::Symbol)
    M, N = size(A)
    K = min(M, N)
    scale, rowdims, coldims = dimsplat(A)
    if d === :Q
        return nodims(N, N)
    elseif d === :L
        return scale ⊗ rowdims ⊗ nodims(K)
    else
        _throw_factor_dimensions(F, d)
    end
end

function _factor_dimensions(F::Union{QR, QRCompactWY}, A::AbstractAxesDimensions, d::Symbol)
    M, N = size(A)
    K = min(M, N)
    scale, rowdims, coldims = dimsplat(A)
    if d === :Q
        return nodims(M, M)
    elseif d === :R
        return scale ⊗ nodims(K) ⊗ coldims
    else
        _throw_factor_dimensions(F, d)
    end
end

function _factor_dimensions(F::QRPivoted, A::AbstractAxesDimensions, d::Symbol)
    p = F.p
    M, N = size(A)
    K = min(M, N)
    scale, rowdims, coldims = dimsplat(A)
    if d === :Q
        return nodims(M, M)
    elseif d === :R
        return scale ⊗ nodims(K) ⊗ coldims[p]
    elseif d === :P
        return _right_permutation_matrix(A, p)
    else
        _throw_factor_dimensions(F, d)
    end
end



for f in (:lowrankupdate, :lowrankupdate!, :lowrankdowndate, :lowrankdowndate!)
    @eval function $f(F::UnitfulFactorization{<:Any,<:Any,<:Cholesky}, v::AbstractUnitfulVector)
        dimsF = dimensions(F); dimsv = dimensions(v)
        if normdims(dimsF)[1] != normdims(dimsv)[1] || dimscale(dimsF) != dimscale(dimsv)^2
            throw(DimensionMismatch("dimension mismatch in $($f)"))
        end
        return UnitfulFactorization($f(values.((F, v))...), dimsF)
    end
end