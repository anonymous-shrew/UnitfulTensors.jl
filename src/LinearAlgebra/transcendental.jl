for f in (:exp ,  :log ,  :cis ,
          :sin ,  :cos ,  :tan ,  :cot ,  :sec ,  :csc ,
         :asin , :acos , :atan , :acot , :asec , :acsc ,
          :sinh,  :cosh,  :tanh,  :coth,  :sech,  :csch,
         :asinh, :acosh, :atanh, :acoth, :asech, :acsch,
          :sind,  :cosd,  :tand,# :coth,  :sech,  :csch, (they don't yet accept matrices in stdlib)
         :asind, :acosd, :atand, :acotd, :asecd, :acscd,
)
    @eval @inline function $f(args::AbstractUnitfulScalarOrTensor...; kwargs...)
        dims = $f(dimensions.(args)...; kwargs...)
        vals = $f(values.(args)...; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return T(vals, dims)
    end
    
    @eval function $f(A::AbstractMatrixDimensions) 
        if !isendomorphic(A)
            throw(DimensionMismatch("$($f) requires an endomorphic matrix"))
        end
        return A
    end
end

for f in (:sincos, :sincosd,
)
    @eval @inline function $f(args::AbstractUnitfulScalarOrTensor...; kwargs...)
        dims = $f(dimensions.(args)...; kwargs...)
        vals = $f(values.(args)...; kwargs...)
        T = promote_unitful($f, vals, dims, args...)
        return T.(vals, dims)
    end
    @eval function $f(A::AbstractMatrixDimensions)
        isendomorphic(A) || throw(DimensionMismatch("$($f) requires an endomorphic matrix"))
        return (A, A)
    end
end

for T in (Number, Irrational{:ℯ})
    @eval @inline function (^)(b::$T, A::AbstractUnitfulScalarOrTensor)
        dims = b^dimensions(A)
        vals = b^values(A)
        T = promote_unitful(^, vals, dims, A)
        return T(vals, dims)
    end
    @eval function (^)(b::$T, A::AbstractMatrixDimensions)
        if !isendomorphic(A)
            throw(DimensionMismatch("only an endomorphic matrix can be an exponent"))
        end
        return A
    end
    @eval function (^)(b::$T, A::AbstractDimensions)
        if A != one(A)
            throw(DimensionMismatch("only a dimensionless quantity can be an exponent"))
        end
        return A
    end
end