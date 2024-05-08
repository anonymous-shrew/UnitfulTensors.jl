# To silence warnings emitted in \ and pinv of inhomogeneous rectangular matrices
loglevel = Base.CoreLogging._min_enabled_level[]
Logging.disable_logging(Logging.Warn)

@testset "addition" begin
    test_unitful_unitless("scalar + scalar", +, ħ, h)
    test_unitful_unitless("vector + vector", +, u1, v1)
    test_unitful_unitless("matrix + matrix", +, A13, B13)
end
@testset "subtraction" begin
    test_unitful_unitless("scalar - scalar", -, ħ, h)
    test_unitful_unitless("vector - vector", -, u1, v1)
    test_unitful_unitless("matrix - matrix", -, A13, B13)
end 
@testset "multiplication" begin
    test_unitful_unitless("scalar * scalar", *, ħ, e)
    test_unitful_unitless("scalar * vector", *, ħ, v1)
    test_unitful_unitless("vector * scalar", *, v1, ħ)
    test_unitful_unitless("scalar * matrix", *, ħ, A13)
    test_unitful_unitless("matrix * scalar", *, A13, ħ)
    test_unitful_unitless("matrix * vector", *, A13, v3)
    test_unitful_unitless("matrix * matrix", *, A13, A32)
    test_unitful_unitless("vector * single-row matrix", *, v3, As1singlerow)
    
    test_unitful_unitless("adjoint * scalar", *, v1', ħ)
    test_unitful_unitless("scalar * adjoint", *, ħ, v1')
    test_unitful_unitless("adjoint * vector", *, v1inv', v1)
    test_unitful_unitless("vector * adjoint", *, v3, v1inv')
    test_unitful_unitless("adjoint * matrix", *, v1inv', A13)
    test_unitful_unitless("single-column matrix * adjoint", *, A3ssinglecol, v1inv')
end
@testset "left division" begin
    test_unitful_unitless("scalar \\ scalar", \, ħ, e)
    test_unitful_unitless("scalar \\ vector", \, ħ, v1)
    test_unitful_unitless("scalar \\ matrix", \, ħ, A13)
    test_unitful_unitless("matrix \\ vector", \, A12square, v1)
    test_unitful_unitless("matrix \\ matrix", \, A12square, A13)
    @testset "overdetermined" begin
        test_unitful_unitless("matrix \\ vector", \, A02tall, v0)
        test_unitful_unitless("matrix \\ matrix", \, A02tall, A03)
        test_unitful_unitless("vector \\ vector", \, u0, v0)
        test_unitful_unitless("vector \\ matrix", \, v0, A03)
    end
    @testset "underdetermined" begin
        test_unitful_unitless("matrix \\ vector", \, A10wide, v1)
        test_unitful_unitless("matrix \\ matrix", \, A10wide, A13)
    end
end
@testset "right division" begin
    test_unitful_unitless("scalar / scalar", /, ħ, e)
    test_unitful_unitless("vector / scalar", /, v1, ħ)
    test_unitful_unitless("matrix / scalar", /, A13, ħ)
    # The next one is probably not useful enough to be implemented
    test_unitful_unitless("vector / vector", /, v1, v3; skip=true)
    test_unitful_unitless("matrix / matrix", /, A32, A12square)
    test_unitful_unitless("vector / single-column matrix", /, v1, A3ssinglecol; skip=true)
    test_unitful_unitless("single-column matrix / vector", /, A3ssinglecol, v1; skip=true)
    test_unitful_unitless("adjoint / matrix", /,  v2inv', A12square)
end 
@testset "inv" begin
    test_unitful("scalar", inv, ħ)
    test_unitful("matrix", inv, A12square)
end
@testset "pinv" begin
    test_unitful("scalar", pinv, ħ)
    test_unitful("vector", pinv, v0)
    test_unitful("matrix", pinv, A12square)
    test_unitful("overdetermined", pinv, A02tall)
    test_unitful("underdetermined", pinv, A10wide)
end
@testset "lmul!" begin
    test_unitful("scalar * matrix", lmul!, ħ, A13)
    test_unitful("matrix * matrix", lmul!, A21upper, A13)
end
@testset "rmul!" begin
    test_unitful("matrix * scalar", rmul!, A13, ħ)
    test_unitful("matrix * matrix", rmul!, A32, A21upper)
end
@testset "ldiv!" begin
    test_unitful("scalar \\ matrix", ldiv!, ħ, A13)
    test_unitful("matrix \\ matrix", ldiv!, A12upper, A13)
    test_unitful("matrix = matrix \\ matrix", ldiv!, A12ldiv13, A12upper, A13)
end
@testset "rdiv!" begin
    test_unitful("matrix / scalar", rdiv!, A13, ħ)
    test_unitful("matrix / matrix", rdiv!, A32, A12upper)
end
@testset "dot" begin
    test_unitful_unitless("scalar ⋅ scalar", ⋅, ħ, e)
    test_unitful_unitless("vector ⋅ vector", ⋅, v1, v1inv)
    test_unitful_unitless("matrix ⋅ matrix", ⋅, A13, A1inv3inv)
    test_unitful("dot(vector, matrix, vector)", ⋅, v1inv, A13, v3)
end
@testset "cross" begin
    test_unitful_unitless("vector × vector", ×, u3d, v3d)
end
@testset "power" begin
    for p in (0, 1, -1, 2, -2, 13, -13, 1//3, -1//3, 1.618, -1.618)
        @testset "x^$p" begin
            kwargs = p isa Integer ? () : (;rtol = 1e-6)
            test_unitful("scalar", x -> x^p, ħ; kwargs...)
            test_unitful("matrix", x -> x^p, A11squareable; kwargs...)
        end
    end
end
@testset "sqrt" begin
    test_unitful("scalar", sqrt, ħ)
    test_unitful("matrix", sqrt, A11squareable)
end
@testset "kron" begin
    for (arg1, desc1) in ((ħ, "scalar"), (v1, "vector"), (A13, "matrix"),
                          (v1', "adjoint"), (A13', "adjoint matrix"))
        for (arg2, desc2) in ((e, "scalar"), (v2inv, "vector"), (A32, "matrix"),
                              (v2inv', "adjoint"), (A32', "adjoint matrix"))
            test_unitful_unitless("kron($desc1, $desc2)", kron, arg1, arg2)
            if !(desc1 == desc2 == "scalar")
                @testset "kron!(..., $desc1, $desc2)" begin
                    dest = kron(arg1, arg2)
                    @test kron!(deepcopy(dest), arg1, arg2) == dest
                end
            end
        end
    end
end

Logging.disable_logging(loglevel - 1)