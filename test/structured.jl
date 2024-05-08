for f in (diag, tril, triu, tril!, triu!)
    @testset "$f" begin
        test_unitful("$f(A)", f, A13)
        test_unitful("$f(A, 1)", x -> f(x, 1), A13)
        test_unitful("$f(A, -1)", x -> f(x, -1), A13)
    end
end
for f in (hermitianpart, hermitianpart!)
    @testset "$f" begin
        test_unitful("$f(A)", f, A11invsym)
        test_unitful("$f(A, :U)", x -> f(x, :U), A11invsym)
        test_unitful("$f(A, :L)", x -> f(x, :L), A11invsym)
    end
end
for f in (Diagonal,)
    @testset "$f" begin
        test_unitful("vector", f, v1)
        test_unitful("matrix", f, A13)
    end
end