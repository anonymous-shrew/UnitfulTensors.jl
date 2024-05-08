@testset "det" begin
    test_unitful("scalar", det, ħ)
    test_unitful("matrix", det, A12square)
end
@testset "logdet" begin
    test_unitful("scalar", logdet, scalar_unitless)
    test_unitful("matrix", logdet, A11endomorphic * A11endomorphic)
end
@testset "logabsdet" begin
    test_unitful("scalar", logabsdet, scalar_unitless)
    test_unitful("matrix", logabsdet, A11endomorphic)
end
@testset "tr" begin
    test_unitful("scalar", tr, ħ)
    test_unitful("matrix", tr, A11squareable)
end
@testset "lyap" begin
    test_unitful("scalar", lyap, ħ, e)
    test_unitful("matrix", lyap, A11squareable, A11invsym)
end
@testset "sylvester" begin
    test_unitful("scalar", sylvester, ħ, h, e)
    test_unitful("matrix", sylvester, A11squareable, A33squareable, A13)
end
for f in (rotate!, reflect!)
    @testset "$f" begin
        test_unitful("$f(..., number, number)", f, v1, u1, real_number, real_number2)
        test_unitful("$f(..., scalar, scalar)", f, v1, u1, real_unitless, real_unitless2)
    end
end
for f in (givens, )
    @testset "$f" begin
        test_unitful("$f(scalar, scalar, ...)", f, real_unitless, real_unitless2, 2, 3)
        test_unitful("$f(vector, ...)", f, v0, 2, 3)
        test_unitful("$f(matrix, ...)", f, A02tall, 2, 3, 3)
    end
end