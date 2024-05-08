prop_test(::typeof(lu), F, A) = F.L * F.U ≈ A[F.p, :]
iter_test(::typeof(lu), F, A) = ((L, U, p) = F; L * U ≈ A[p, :])
prop_test(::typeof(bunchkaufman), F, A) = inv(F.P) * F.U * F.D * F.U' * inv(F.P') ≈ A
iter_test(::typeof(bunchkaufman), F, A) = ((D, U, p) = F; U * D * U' ≈ A[p, p])
prop_test(::typeof(cholesky), F, A) = F.L * F.U ≈ A
iter_test(::typeof(cholesky), F, A) = ((L, U) = F; L * U ≈ A)
prop_test(::typeof(ldlt), F, A) = F.L * F.D * F.L' ≈ A
iter_test(::typeof(ldlt), F, A) = true # LDLt is not iterable
prop_test(::typeof(eigen), F, A) = ((vals, vecs) = (F.values, F.vectors); vecs * Diagonal(vals) ≈ A * vecs)
iter_test(::typeof(eigen), F, A) = ((vals, vecs) = F; vecs * Diagonal(vals) ≈ A * vecs)
prop_test(::typeof(hessenberg), F, A) = F.Q * F.H * inv(F.Q) ≈ A
iter_test(::typeof(hessenberg), F, A) = ((Q, H, μ) = F; Q * H * inv(Q) ≈ A)
prop_test(::typeof(schur), F, A) = F.vectors * F.Schur * inv(F.vectors) ≈ A
iter_test(::typeof(schur), F, A) = ((T, Z, vals) = F; Z * T * inv(Z) ≈ A)
prop_test(::typeof(svd), F, A) = F.U * Diagonal(F.S) * F.Vt ≈ A
iter_test(::typeof(svd), F, A) = ((U, S, V) = F; U * Diagonal(S) * V' ≈ A)
prop_test(::typeof(lq), F, A) = F.L * F.Q[1:size(F.L, 2), :] ≈ A
iter_test(::typeof(lq), F, A) = ((L, Q) = F; L * Q[1:size(L, 2), :] ≈ A)
prop_test(::typeof(qr), F, A) = F.Q[:, 1:size(F.R, 1)] * F.R ≈ A
iter_test(::typeof(qr), F, A) = ((Q, R) = F; Q[:, 1:size(R, 1)] * R ≈ A)

# Generalized factorizations
prop_test(::typeof(eigen), F, A, B) = ((vals, vecs) = (F.values, F.vectors); B * vecs * Diagonal(vals) ≈ A * vecs)
iter_test(::typeof(eigen), F, A, B) = ((vals, vecs) = F; B * vecs * Diagonal(vals) ≈ A * vecs)
prop_test(::typeof(schur), F, A, B) = F.Q * F.S * F.Z' ≈ A && F.Q * F.T * F.Z' ≈ B &&
                                      all(rank(values(βi*A - αi*B), # testing that α ./ β are generalized eigenvalues
                                               atol = 1e-15 * maximum((A, B) .|> values .|> norm)
                                              ) < size(A, 1) for (αi, βi) in zip(F.α, F.β))
iter_test(::typeof(schur), F, A, B) = ((S, T, Q, Z, α, β) = F; Q * S * Z' ≈ A && Q * T * Z' ≈ B) &&
                                       all(rank(values(βi*A - αi*B), # testing that α ./ β are generalized eigenvalues
                                                atol = 1e-15 * maximum((A, B) .|> values .|> norm)
                                               ) < size(A, 1) for (αi, βi) in zip(α, β))
prop_test(::typeof(svd), F, A, B) = F.U * F.D1 * F.R0 * F.Q' ≈ A && F.V * F.D2 * F.R0 * F.Q' ≈ B
iter_test(::typeof(svd), F, A, B) = ((U, V, Q, D1, D2, R0) = F; U * D1 * R0 * Q' ≈ A && V * D2 * R0 * Q' ≈ B)

function _sort_eigvals(vals)
    s = sort(vals, by = reim ∘ value)
    return vals isa UnitfulTensor ? UnitfulTensor(s) : s
end

test_unitful_factorization(lu, A13)
@testset "cholesky and friends" begin
    test_unitful_factorization(bunchkaufman, -A11invposdef)
    test_unitful_factorization(cholesky, A11invposdef)
    test_unitful_factorization(ldlt, A11invsymtridiag)
    for (f, s) in zip((lowrankupdate, lowrankupdate!, lowrankdowndate, lowrankdowndate!), (1, 1, -1, -1))
        @testset "$f" begin
            @test all((f(cholesky(A11invposdef), deepcopy(v1sqrt))..., )
                       .≈ (cholesky(A11invposdef + s*v1sqrt*v1sqrt')..., ))
        end
    end
end
@testset "eigen and friends" begin
    test_unitful_factorization(eigen, A11squareable)
    test_unitful_factorization(hessenberg, A11squareable)
    test_unitful_factorization(schur, A11squareable)
    test_unitful_factorization(svd, A04homogeneous)
    test_unitful("eigvals", eigvals, A11squareable)
    test_unitful("eigvals!", eigvals!, A11squareable)
    @test equal_sametype(eigvecs(A11squareable), eigen(A11squareable).vectors)
    test_unitful("eigmax", eigmax, A11realeigvals)
    test_unitful("eigmin", eigmin, A11realeigvals)
    test_unitful("svdvals", svdvals, A04homogeneous)
    test_unitful("svdvals!", svdvals!, A04homogeneous)
    @testset "generalized" begin
        test_unitful_factorization(eigen, A11squareable, B11squareable)
        test_unitful_factorization(schur, A11squareable, B11squareable)
        test_unitful_factorization(svd, A04homogeneous, A44squarehomogeneous)
        test_unitful("eigvals", _sort_eigvals ∘ eigvals, A11squareable, B11squareable)
        test_unitful("eigvals!", _sort_eigvals ∘ eigvals!, A11squareable, B11squareable)
        @test equal_sametype(eigvecs(A11squareable, B11squareable), eigen(A11squareable, B11squareable).vectors)
        test_unitful("svdvals", svdvals, A04homogeneous, A44squarehomogeneous)
        test_unitful("svdvals!", svdvals!, A04homogeneous, A44squarehomogeneous)
    end
end
@testset "qr and friends" begin
    test_unitful_factorization(lq, A10wide)
    test_unitful_factorization(qr, A02tall)
end