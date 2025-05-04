using TensorOperations

@testset "TensorOperations" begin
    test_unitful("addition", (A, B) -> (@tensoropt Z[i, j] := A[i, j] + 2 * B[i, j]; Z), A13, B13)
    test_unitful("matrix multiplication", (A, C) -> (@tensoropt Z[i, k] := A[i, j] * C[j, k]; Z), A13, A32)
    test_unitful("general contraction", (A, B, C, D) -> (@tensoropt Z[i, j, k] := A[i, i'] * B[j, j'] * C[k, k'] * D[i', j', k']; Z), A11squareable, A32, A13, A123)
    test_unitful("trace", A -> (@tensoropt Z = A[i, i]; Z), A11squareable)
    test_unitful("outer product & index permutation", (A, C) -> (@tensoropt Z[i, j, k, l] := A[j, l] * C[k, i]; Z), A13, A32)
end