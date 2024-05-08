for f in (exp ,  log ,
          sin ,  cos ,  tan ,  cot ,  sec ,  csc ,
         asin , acos , atan , acot , asec , acsc ,
          sinh,  cosh,  tanh,  coth,  sech,  csch,
         asinh, acosh, atanh, acoth, asech, acsch,
          sind,  cosd,  tand,# coth,  sech,  csch, (they don't yet accept matrices in Base)
         asind, acosd, atand, acotd, asecd, acscd,
         sincos, sincosd
)
    test_unitful("$f", f, A11endomorphic)
end

@testset "exponentiation" begin
    for b in (1, 2, 1//3, 1.618, 1.618im, ℯ)
        @testset "$b^x" begin
            test_unitful("scalar", x -> b^x, h/ħ)
            test_unitful("matrix", x -> b^x, A11endomorphic)
        end
    end
end