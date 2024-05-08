for f in (:det, :logdet, :tr, :lyap, :sylvester)
    @eval @inline function $f(args::AbstractUnitfulScalarOrTensor...; kwargs...)
        dims = $f(dimensions.(args)...; kwargs...)
        vals = $f(values.(args)...; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return T(vals, dims)
    end
end

for f in (:logabsdet, )
    @eval @inline function $f(args::AbstractUnitfulScalarOrTensor...; kwargs...)
        dims = $f(dimensions.(args)...; kwargs...)
        vals, sign = $f(values.(args)...; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return T(vals, dims), sign
    end
end

function det(A::AbstractMatrixDimensions)
    N = checksquare(A)
    scale, rowdims, coldims = dimsplat(A)
    return scale^N * prod(rowdims) * prod(coldims)
end

det(A::AbstractDimensions) = A

logdet(A::AbstractDimensionsOrAxesDimensions) = log(det(A))
logabsdet(A::AbstractDimensionsOrAxesDimensions) = log(det(A))

function tr(A::AbstractMatrixDimensions) 
    if !issquareable(A)
        throw(DimensionMismatch("tr requires a squareable matrix"))
    else
        return dimscale(A)
    end
end

function sylvester(A::AbstractMatrixDimensions, B::AbstractMatrixDimensions, C::AbstractMatrixDimensions)
    if dimscale(A) != dimscale(B) ||
        !issquareable(A) || !issquareable(B) ||
        normdims(A)[1] != normdims(C)[1] || normdims(B)[2] != normdims(C)[2]
            throw(DimensionMismatch("incompatible physical dimensions in sylvester"))
    else        
        return C / dimscale(A)
    end
end

sylvester(a::AbstractDimensions, b::AbstractDimensions, c::AbstractDimensions) = -c / (a + b)

function lyap(A::AbstractMatrixDimensions, C::AbstractMatrixDimensions)
    if  !issquareable(A) || !issymmetric(C) || normdims(A)[1] != normdims(C)[1]
        throw(DimensionMismatch("incompatible physical dimensions in lyap"))
    else        
        return C / dimscale(A)
    end
end

lyap(a::AbstractDimensions, c::AbstractDimensions) = c / a

for f in (:rotate!, :reflect!)
    @eval @inline function $f(x::AbstractUnitfulVector, y::AbstractUnitfulVector, c, s)
        args = (x, y, c, s); kwargs = ()
        dims = $f(dimensions.(args)...; kwargs...)
        vals = $f(values.(args)...; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return T.(vals, dims)
    end
    @eval function $f(x::AbstractVectorDimensions, y::AbstractVectorDimensions, c, s)
        if x != y || c != one(c) || s != one(s)
            throw(DimensionMismatch("$($f) requires vectors with the same physical dimensions\
                                     and dimensionless c, s"))
        end
        return x, y
    end
end

for f in (:givens, )
    @eval @inline function $f(f::AbstractUnitfulScalar, g::AbstractUnitfulScalar, i1::Integer, i2::Integer)
        args = (f, g, i1, i2); kwargs = ()
        dimsf, dimsg = dimensions.((f, g))
        if dimsf != dimsg
            throw(DimensionMismatch("$($f) requires scalars with the same physical dimensions"))
        end
        dims = dimsf
        G, vals = $f(values.(args)...; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return G, T(vals, dims)
    end
end