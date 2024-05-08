atol = 0; rtol = 1

@testset "cond" begin
    test_unitful("cond(A)", cond, A04homogeneous)
    test_unitful("cond(A, Inf)", x -> cond(x, Inf), A44squarehomogeneous)
end
@testset "condskeel" begin
    test_unitful("condskeel(A)", condskeel, A14squarehomogeneousrows)
    test_unitful("condskeel(A, 1)", x -> condskeel(x, 1), A14squarehomogeneousrows)
    test_unitful("condskeel(A, v)", condskeel, A14squarehomogeneousrows, v4)
    test_unitful("condskeel(A, v, 2)", (x, y) -> condskeel(x, y, 2), A14squarehomogeneousrows, v4)
end
for f in (norm, normalize, normalize!)
    @testset "$f" begin
        if f !== normalize!
            test_unitful("scalar", f, ħ)
            test_unitful("scalar, p = -0.5", x -> f(x, -0.5), ħ)
        end
        test_unitful("vector", f, v0)
        test_unitful("vector, p = -0.5", x -> f(x, -0.5), v0)
        test_unitful("matrix", f, A04homogeneous)
        test_unitful("matrix, p = -0.5", x -> f(x, -0.5), A04homogeneous)
    end
end
@testset "opnorm" begin
    test_unitful("scalar", opnorm, ħ)
    test_unitful("matrix", opnorm, A04homogeneous)
    test_unitful("scalar, p = 1", x -> opnorm(x, 1), ħ)
    test_unitful("matrix, p = 1", x -> opnorm(x, 1), A04homogeneous)
end
@testset "rank" begin
    @testset "rank(A)" begin
        @test rank(A13) == rank(values(A13))
    end
    @testset "rank(A; atol, rtol)" begin
        @test rank(A13; atol = atol, rtol = rtol) == rank(values(A13); atol = atol, rtol = rtol)
    end
end
@testset "nullspace" begin
    @testset "nullspace(A)" begin
        @test nullspace(A13) isa UnitfulScalarOrTensor
        @test values(nullspace(A13)) == nullspace(values(A13))
        @test (A13 * nullspace(A13); true) # testing multipliability
    end
    @testset "nullspace(A; atol, rtol)" begin
        @test nullspace(A13) isa UnitfulScalarOrTensor
        @test values(nullspace(A13; atol = atol, rtol = rtol)) == nullspace(values(A13); atol = atol, rtol = rtol)
        @test (A13 * nullspace(A13; atol = atol, rtol = rtol); true) # testing multipliability
    end
end