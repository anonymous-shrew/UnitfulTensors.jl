unitfulize(x) = convert(Unitful.Dimensions, x)
double_convert(x) = SIDimensions(unitfulize(x))

rand_num = randn()
rand_dims = prod((𝐋, 𝐌, 𝐓, 𝐈, 𝚯, 𝐍, 𝐉) .^ randn(7))
unitful_rand_dims = unitfulize(rand_dims)

@testset "scalar interface" begin
    @test one(𝐋) === NoDims
    @test value(rand_num) === rand_num
    @test values(rand_num) === rand_num
    @test dimensions(rand_num) === NoDims
    @test value(ħ) === valħ
    @test dimensions(ħ) === 𝐋^2 * 𝐌 * 𝐓^-1
    @test dimensions(2 * u"m/m") === NoDims
end

@testset "scalar conversion" begin
    @test rand_dims === double_convert(rand_dims)
    @test SIDimensions(unitful_rand_dims) === convert(SIDimensions, unitful_rand_dims)
end

@testset "scalar show" begin
    @test repr(rand_dims) === repr(convert(Unitful.Dimensions, rand_dims))
    @test contains(repr(ħ), repr(value(ħ)))
    @test contains(repr(ħ), repr(Unitful.upreferred(unitfulize(dimensions(ħ)))))
    @test [contains(repr(ħ), x) for x in ("kg", "m", "s")] == fill(true, 3)
end

@testset "arithmetic with AbstractDimensions" begin
    @test 𝐋^1000000 / 𝐋^999999.5 === sqrt(𝐋)
    @test sqrt(𝐋)^2000000 * sqrt(𝐋)^(-1999999) === sqrt(𝐋)
    @test sqrt(𝐋)^2000000 * sqrt(𝐋)^(-1999998) === 𝐋
    @test rand_dims^1000 / rand_dims^999 / rand_dims ≈ NoDims
    @test rand_dims^1000 ≈ rand_dims^500 * rand_dims^500
    @test rand_dims^1000 == rand_dims^1000. == rand_dims^(1000//1)
    @test inv(rand_dims) == NoDims / rand_dims
end