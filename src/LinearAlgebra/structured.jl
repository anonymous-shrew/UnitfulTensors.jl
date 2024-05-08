for f in (:diag, )
    @eval @inline function $f(A::AbstractUnitfulScalarOrTensor, k::Integer=0)
        dims = $f(dimensions(A), k)
        vals = $f(values(A), k)
        T = promote_unitful($f, vals, dims, A)
        return T(vals, dims)
    end
end

# This is redundant since Julia 1.11, when diagind(A, k, IndexCartesian()) was introduced
function diag(A::AbstractMatrixDimensions, k::Integer=0)
    scaleA, rowdims, coldims = dimsplat(A)
    m, n = size(A)
    rowrange = (1:m) ∩ ((1:n) .- k)
    colrange = ((1:m) .+ k) ∩ (1:n)
    if k == 0
        scale = scaleA
        @views dims = (rowdims[rowrange] .* coldims[colrange], )
    else
        s = (k > 0) ? coldims[1 + k] : rowdims[1 - k]
        scale = scaleA * s
        @views dims = (rowdims[rowrange] .* coldims[colrange] ./ s, )
    end
    T = promote_dims(diag, (dims, scale), A)
    return T(dims, scale; normalize = false)
end

for f in (:tril, :tril!, :triu, :triu!)
    @eval function $f(A::AbstractUnitfulMatrix)
        vals = $f(values(A)); dims = dimensions(A)
        T = promote_unitful($f, vals, dims, A)
        return T(vals, dims)
    end
    @eval function $f(A::AbstractUnitfulMatrix, k::Integer)
        vals = $f(values(A), k); dims = dimensions(A)
        T = promote_unitful($f, vals, dims, A)
        return T(vals, dims)
    end
end

for f in (:hermitianpart, :hermitianpart!)
    @eval @inline function $f(A::AbstractUnitfulScalarOrTensor; kwargs...)
        args = (A, )
        dims = $f(dimensions.(args)...; kwargs...)
        vals = $f(values.(args)...; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return T(vals, dims)
    end
    @eval @inline function $f(A::AbstractUnitfulScalarOrTensor, optarg::Symbol; kwargs...)
        args = (A, )
        dims = $f(dimensions.(args)..., optarg; kwargs...)
        vals = $f(values.(args)..., optarg; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return T(vals, dims)
    end
    @eval function $f(A::AbstractMatrixDimensions, uplo::Symbol=:U)
        if !issymmetric(A)
            throw(DimensionMismatch("$($f) requires a dimensionally symmetric matrix"))
        end
        return A
    end
    @eval $f(A::AbstractDimensions, maybe_optarg...) = A
end

function Diagonal(A::AbstractUnitfulTensor)
    dims = Diagonal(dimensions(A))
    vals = Diagonal(values(A))
    T = promote_unitful(Diagonal, vals, dims, A)
    return T(vals, dims)
end

function Diagonal(A::AbstractVectorDimensions)
    scale = dimscale(A)
    dims = AxisDimensions(sqrt.(normdims(A)[1]))
    T = promote_dims(Diagonal, ((dims, dims), scale), A)
    return T((dims, dims), scale)
end

function Diagonal(A::AbstractMatrixDimensions)
    m, n = size(A)
    if m == n
        return A
    elseif m > n
        return A[1:n, :]
    else
        return A[:, 1:m] 
    end
end