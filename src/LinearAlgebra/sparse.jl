@inline function sparse(A::AbstractUnitfulMatrix)
    dims = sparse(dimensions(A))
    vals = sparse(values(A))
    T = promote_unitful(sparse, vals, dims, A)
    return T(vals, dims)
end

sparse(A::AbstractAxesDimensions) = A