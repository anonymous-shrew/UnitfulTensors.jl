dimħ = dimensions(ħ)
Z = randomUnitfulTensor((), dimħ)
V = randomUnitfulTensor((dims0, ), dimħ)
M = randomUnitfulTensor((dims0, dims0), dimħ)

testobjects = (("scalar", ħ), ("0d", Z), ("vector", V), ("matrix", M), ("adjoint", V'))

@testset "broadcasting" begin
    @testset "addition" begin
        for (desc1, x1) in testobjects
            for (desc2, x2) in testobjects
                test_unitful_unitless("$desc1 .+ $desc2", (x, y) -> x .+ y, x1, x2)
            end
        end
    end
end