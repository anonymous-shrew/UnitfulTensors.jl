basedims = [𝐋, 𝐌]
dims = ([NoDims, 𝐋, 𝐌^(-1), 𝐋/𝐌], ([NoDims, dim, dim^2] for dim in [basedims; basedims.^(-1)])...)
scale = 𝚯
U = randomUnitfulTensor(dims, scale)
U2 = randomUnitfulTensor(dims, scale)
A = dimensions(U)

I0 = (1, 2, 3, 1, 2)
I1 = ([1 2; 3 4], fill(2), 3, BigInt(2), 3:-2:1)
I2 = (CartesianIndex(1, 3),
      CartesianIndex(), :, CartesianIndices(()), :, fill(CartesianIndex(), 2),
      CartesianIndex(2))
I3 = (CartesianIndices((1:3, 3:-2:1)), CartesianIndices((2:2, )),
     [CartesianIndex(1, 1) CartesianIndex(2, 2)
      CartesianIndex(2, 2) CartesianIndex(3, 3)],
      [1], 1, true:true)
I4 = (CartesianIndices((1:3, 3:-2:1)),
      [CartesianIndex(1, 1) CartesianIndex(2, 2)
       CartesianIndex(2, 2) CartesianIndex(3, 3)],
      [false, true, true], [1], 1, true:true)
I5 = ([], :, CartesianIndices((2:1, 3:1)), [;;])
I6 = (rand(Bool, size(A)), )

array(A) = Array(A)
array(A::AbstractDimensions) = A
array(A::AbstractUnitfulScalar) = A

@testset "getindex" begin
    @testset "AxesDimensions" for I in (I0, I1, I2, I3, I4, I5, I6)
        @test array(A[I...]) == array(A)[I...]
        @test A[I...] isa Union{AxesDimensions, SIDimensions}
    end
    @testset "UnitfulTensor" for I in (I0, I1, I2, I3, I4, I5, I6)
        @test array(U[I...]) == array(U)[I...]
        @test U[I...] isa Union{UnitfulTensor, UnitfulScalar}
    end
end

Nd_to_1d(x) = x
Nd_to_1d(x::AbstractArray) = vec(x)
Nd_to_1d(x::CartesianIndices) = x
Nd_to_1d(x::AbstractArray{Bool}) = x

@testset "setindex!" begin
    @testset "AxesDimensions" for I in (I0, I1, I2, I3, I4, I5, I6)
        # because A1[I...] = A1[I...] will fail with multidimensional arrays as indices
        # this is a bug (or a feature?) of Julia Base
        J = Nd_to_1d.(I)
        @test (A1 = deepcopy(A); A1[J...] = A1[J...]; A == A1)
    end
    @testset "UnitfulTensor" for I in (I0, I1, I2, I3, I4, I5, I6)
        J = Nd_to_1d.(I)
        @test (U1 = deepcopy(U); U1[J...] = U2[J...]; dimensions(U1) == dimensions(U))
        @test (U1 = deepcopy(U); U1[J...] = U2[J...]; values(U1)[J...] == values(U2)[J...])
    end
end