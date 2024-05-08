using UnitfulTensors, LinearAlgebra
using UnitfulTensors: promote_unitful
using LinearAlgebra: AdjOrTrans, Givens
import Unitful, Logging
using Test
using Random: seed!
import Base: values

seed!(0)



function values(x::AbstractUnitfulScalar, conversion_factors::NTuple{7, Number})
    val = values(x)
    dims = dimensions(x)
    return val * prod(conversion_factors .^ dimexps(dims))
end

function values(x::AbstractUnitfulTensor, conversion_factors::NTuple{7, Number})
    dest = _similar(values(x))
    return copyto!(dest, values.(x, (conversion_factors, )))
end

_similar(x) = similar(x)
_similar(x::Hermitian) = similar(parent(x))

values(x, conversion_factors::NTuple{7, Number}) = values(x)



equal_sametype(x, y) = typeof(x) === typeof(y) && x == y

_dims_type(x) = AbstractDimensions
_dims_type(x::AbstractArray) = AbstractAxesDimensions{ndims(x)}
_dims_type(x::AdjOrTrans) = AdjointAxesDimensions{<:Any,
                            <:AbstractAxesDimensions{ndims(parent(x))}}

_isapprox(x, y; kwargs...) = ≈(x, y; kwargs...)
_isapprox(x::Givens, y::Givens; kwargs...) = (x.i1, x.i2) == (y.i1, y.i2) && all(
                                             (≈).((x.c, x.s), (y.c, y.s); kwargs...))

ξ = (abs.(randn(7))..., )

function test_unitful_nonmutating_single_return(desc, f, args...; skip=false, kwargs...)
    @testset "$desc" begin
        @test f(args...) isa UnitfulScalarOrTensor || equal_sametype(f(args...), f(values.(args)...))
        @test values(f(args...)) == f(values.(args)...)
        @test typeof(values(f(args...))) === typeof(f(values.(args)...))
        if f(args...) isa UnitfulScalarOrTensor
            @test typeof(dimensions(f(args...))) <: _dims_type(f(values.(args)...))
        end
        @test _isapprox(values(f(args...), ξ), f(values.(args, (ξ, ))...); kwargs...)
    end
end

function test_unitful(desc, f, args...; skip=false, kwargs...)
    g(xs...) = f(deepcopy.(xs)...)
    @testset "$desc" begin
        if skip
            @test true skip=true
        else
            res = g(args...)
            if res isa Tuple
                for i in 1:length(res)
                    h(xs...) = g(xs...)[i]
                    test_unitful_nonmutating_single_return("($desc)[$i]", h, args...; kwargs...)
                end
            else
                test_unitful_nonmutating_single_return("$desc", g, args...; kwargs...)
            end
        end
    end
end

function _strip_dims(A)
    vals = values(A)
    dims = dimensions(vals) # nodims(...)
    T = promote_unitful(_strip_dims, vals, dims, A)
    return T(vals, dims)
end

function test_unitful_unitless(desc, f, A, B; kwargs...)
    @testset "$desc" begin
        test_unitful("unitful, unitful", f, A, B; kwargs...)
        test_unitful("unitful, unitless", f, _strip_dims(A), values(B); kwargs...)
        test_unitful("unitless, unitful", f, values(A), _strip_dims(B); kwargs...)
    end
end

function test_unitful_factorization(desc::AbstractString, f, prop_test, iter_test, args...; skip=false)
    g(xs...) = f(deepcopy.(xs)...)
    factors(xs...) = (f ∉ (ldlt, ldlt!)) ? (g(xs...)..., ) : (g(xs...).L, g(xs...).D) # LDLt is not iterable
    @testset "$desc" begin
        @test g(args...) isa Union{UnitfulFactorization, UnitfulGeneralizedFactorization} skip=skip
        @test all(u isa UnitfulScalarOrTensor || u == v
            for (u, v) in zip(factors(args...), factors(values.(args)...))) skip=skip
        @test values.(factors(args...)) == factors(values.(args)...) skip=skip
        @test prop_test(g(args...), args...) skip=skip
        @test iter_test(g(args...), args...) skip=skip
    end
end

function test_unitful_factorization(f, args...; skip=false)
    p = (F, As...) -> prop_test(f, F, As...)
    i = (F, As...) -> iter_test(f, F, As...)
    f! = eval(Symbol(f, !))
    @testset "$f?!" begin
        test_unitful_factorization("$f", f, p, i, args...; skip=skip)
        test_unitful_factorization("$f!", f!, p, i, args...; skip=skip)
    end
end



randomdimensions() = SIDimensions(Tuple(rand(-10:10, 7) .// rand(2 .^ (0:4), 7)))

function randomUnitfulTensor(dims, scale::AbstractDimensions)
    axesdims = AxesDimensions(dims, scale)
    vals = randn(size(axesdims))
    return UnitfulTensor(vals, axesdims)
end

randomUnitfulTensor(dims) = randomUnitfulTensor(dims, randomdimensions())



valħ = 1.054571817e-34
ħ = valħ * u"J*s"
h = 6.62607015e-34 * u"J*s"
e = 1.602176634e-19 * u"C"
real_number = randn()
real_number2 = randn()
real_unitless = real_number * u"s/s"
real_unitless2 = real_number2 * u"s/s"
scalar_unitless = randn(Complex{Float64}) * u"s/s"

dimss = [NoDims]
dims0 = fill(NoDims, 13)
dims1 = [NoDims, 𝐌^(-1)/𝐋, 𝐓^(1/2)/𝐋]
dims2 = [NoDims, 𝚯*𝐍^2/𝐉^3.5, 𝐋/𝐓*𝐈^(-1/4)]
dims3 = [NoDims, 𝐓, 𝐋^(-1), 𝐌^(1/2)]
dims3d = fill(NoDims, 3)
dims4 = fill(NoDims, length(dims1))

A12square = randomUnitfulTensor((dims1, inv.(dims2)))
A13 = randomUnitfulTensor((dims1, inv.(dims3)))
B13 = randomUnitfulTensor((dims1, inv.(dims3)), dimscale(A13))
A32 = randomUnitfulTensor((dims3, inv.(dims2)))
A1inv3inv = randomUnitfulTensor((inv.(dims1), dims3))
A11squareable = randomUnitfulTensor((dims1, inv.(dims1)))
B11squareable = randomUnitfulTensor((dims1, inv.(dims1)))
A33squareable = randomUnitfulTensor((dims3, inv.(dims3)), dimscale(A11squareable))
A10wide = randomUnitfulTensor((dims1, inv.(dims0)))
A02tall = randomUnitfulTensor((dims0, inv.(dims2)))
A03 = randomUnitfulTensor((dims0, inv.(dims3)))
A04homogeneous = randomUnitfulTensor((dims0, inv.(dims4)))
A44squarehomogeneous = randomUnitfulTensor((dims4, inv.(dims4)), dimscale(A04homogeneous))
A14squarehomogeneousrows = randomUnitfulTensor((dims1, inv.(dims4)))
A11endomorphic = randomUnitfulTensor((dims1, inv.(dims1)), NoDims)
A11invsym = randomUnitfulTensor((dims1, dims1))
A11invposdef = A10wide * A10wide'
A11invsymtridiag = UnitfulTensor(-SymTridiagonal(values(A11invposdef)), dimensions(A11invposdef))
As1singlerow = randomUnitfulTensor((dimss, inv.(dims1)))
A3ssinglecol = randomUnitfulTensor((dims3, inv.(dimss)))
A21square = randomUnitfulTensor((dims2, inv.(dims1)))
A12upper = UnitfulTensor(UpperTriangular(values(A12square)), dimensions(A12square))
A21upper = UnitfulTensor(UpperTriangular(values(A21square)), dimensions(A21square))
A12ldiv13 = randomUnitfulTensor((dims2, inv.(dims3)), dimscale(A12square) \ dimscale(A13))
A11realeigvals = UnitfulTensor((x -> x + x')(values(A11squareable)), dimensions(A11squareable))

v0 = randomUnitfulTensor((dims0, ), (𝐋/𝐓)^10)
u0 = randomUnitfulTensor((dims0, ), (𝐋/𝐓)^10)
v1 = randomUnitfulTensor((dims1, ), (𝐋/𝐓)^10)
u1 = randomUnitfulTensor((dims1, ), (𝐋/𝐓)^10)
v3 = randomUnitfulTensor((dims3, ), (𝐋/𝐓)^10)
v1inv = randomUnitfulTensor((inv.(dims1), ), (𝐋/𝐓)^10)
v2inv = randomUnitfulTensor((inv.(dims2), ), (𝐋/𝐓)^10)
u3d = randomUnitfulTensor((dims3d, ), (𝐋/𝐓)^10)
v3d = randomUnitfulTensor((dims3d, ), (𝐋/𝐓)^10)
v4 = randomUnitfulTensor((dims4, ), (𝐋/𝐓)^10)
v1sqrt = A10wide[:, 1]



function test_files(names...)
    for name in names
        try
            @testset "$name" begin
                include("$name.jl")
            end
        catch
        end
    end
end

@testset "all" begin
    test_files("subtyping")
    test_files("fastquantities")
    test_files("indexing")
    test_files("arithmetic")
    test_files("transcendental")
    test_files("generic")
    test_files("norm")
    test_files("factorizations")
    test_files("structured")
end