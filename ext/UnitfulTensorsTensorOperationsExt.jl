module UnitfulTensorsTensorOperationsExt

import Base: *

using UnitfulTensors
using UnitfulTensors: allmap, promote_unitful, promote_dims

using TensorOperations: TupleTools, VectorInterface,
                        DefaultAllocator, Index2Tuple,
                        linearize,
                        argcheck_tensoradd, argcheck_tensorcontract, argcheck_tensortrace

import TensorOperations: scalartype,
                         dimcheck_tensoradd, dimcheck_tensorcontract, dimcheck_tensortrace,
                         tensoralloc_add, tensoralloc_contract,
                         tensoradd!, tensorcontract!, tensortrace!

# This fixes promote_contract and friends
*(x::AbstractUnitfulScalar, ::VectorInterface.One) = x
*(x::AbstractDimensions, ::VectorInterface.One) = x

VectorInterface.scalartype(::Type{T}) where T <: AbstractUnitfulScalar = T
VectorInterface.scalartype(::Type{T}) where T <: AbstractDimensions = T

######################### Allocation of UnitfulTensors #########################

function tensoralloc_add(TC::Type{<:AbstractUnitfulScalar},
                         A, pA::Index2Tuple, conjA::Bool,
                         istemp::Val=Val(false), allocator=DefaultAllocator())
    dims =                 tensoralloc_add(  dimensions(TC), dimensions(A), pA, conjA, istemp, allocator)
    vals =                 tensoralloc_add(      values(TC),     values(A), pA, conjA, istemp, allocator)
       T = promote_unitful(tensoralloc_add, vals, dims, TC ,            A , pA, conjA, istemp, allocator)
    return T(vals, dims)
end

function tensoralloc_contract(TC::Type{<:AbstractUnitfulScalar},
                              A, pA::Index2Tuple, conjA::Bool,
                              B, pB::Index2Tuple, conjB::Bool,
                              pAB::Index2Tuple,
                              istemp::Val=Val(false), allocator=DefaultAllocator())
    dims =                 tensoralloc_contract(  dimensions(TC), dimensions(A), pA, conjA, dimensions(B), pB, conjB, pAB, istemp, allocator)
    vals =                 tensoralloc_contract(      values(TC),     values(A), pA, conjA,     values(B), pB, conjB, pAB, istemp, allocator)
    T = promote_unitful(tensoralloc_contract, vals, dims, TC ,            A , pA, conjA,            B , pB, conjB, pAB, istemp, allocator)
    return T(vals, dims)
end

######################### Allocation of AxesDimensions #########################

function tensoralloc_add(TC::Type{<:AbstractDimensions},
                         A, pA::Index2Tuple, conjA::Bool,
                         istemp::Val=Val(false), allocator=DefaultAllocator())
    dims = map(n -> normdims(A)[n], linearize(pA))
    scale = dimscale(A)
    T = promote_dims(tensoralloc_add, (dims, scale), TC, A, pA, conjA, istemp, allocator)
    return T(dims, scale)
end

function tensoralloc_contract(TC::Type{<:AbstractDimensions},
                              A, pA::Index2Tuple, conjA::Bool,
                              B, pB::Index2Tuple, conjB::Bool,
                              pAB::Index2Tuple,
                              istemp::Val=Val(false), allocator=DefaultAllocator())
    lA = length(pA[1])
    dims = map(n -> n <= lA ? normdims(A)[pA[1][n]] : normdims(B)[pB[2][n - lA]], linearize(pAB))
    scale = dimscale(A) * dimscale(B)
    T = promote_dims(tensoralloc_contract, (dims, scale), TC, A, pA, conjA, B, pB, conjB, pAB, istemp, allocator)
    return T(dims, scale)
end

######################### Operations on UnitfulTensors #########################

function tensoradd!(C::AbstractUnitfulTensor,
                    A, pA::Index2Tuple, conjA::Bool,
                    α::Number, β::Number,
                    backend, allocator)
    tensoradd!(dimensions(C), dimensions(A), pA, conjA, dimensions(α), dimensions(β), backend, allocator)
    tensoradd!(    values(C),     values(A), pA, conjA,     values(α),     values(β), backend, allocator)
    return C
end

function tensortrace!(C::AbstractUnitfulTensor,
                      A, p::Index2Tuple, q::Index2Tuple, conjA::Bool,
                      α::Number, β::Number,
                      backend, allocator)
    tensortrace!(dimensions(C), dimensions(A), p, q, conjA, dimensions(α), dimensions(β), backend, allocator)
    tensortrace!(    values(C),     values(A), p, q, conjA,     values(α),     values(β), backend, allocator)
    return C
end

function tensorcontract!(C::AbstractUnitfulTensor,
                         A, pA::Index2Tuple, conjA::Bool,
                         B, pB::Index2Tuple, conjB::Bool,
                         pAB::Index2Tuple,
                         α::Number, β::Number,
                         backend, allocator)
    tensorcontract!(dimensions(C), dimensions(A), pA, conjA, dimensions(B), pB, conjB, pAB, dimensions(α), dimensions(β), backend, allocator)
    tensorcontract!(    values(C),     values(A), pA, conjA,     values(B), pB, conjB, pAB,     values(α),     values(β), backend, allocator)
    return C
end

######################### Operations on AxesDimensions #########################

function tensoradd!(C::AbstractAxesDimensions,
                    A, pA::Index2Tuple, conjA::Bool,
                    α::AbstractDimensions, β::AbstractDimensions,
                    backend, allocator)
    β == one(β) || throw(DimensionMismatch("β must be dimensionless"))
    argcheck_tensoradd(C,   A, pA)
    dimcheck_tensoradd(C, α*A, pA)
    return C
end

function tensortrace!(C::AbstractAxesDimensions,
                      A, p::Index2Tuple, q::Index2Tuple, conjA::Bool,
                      α::AbstractDimensions, β::AbstractDimensions,
                      backend, allocator)
    β == one(β) || throw(DimensionMismatch("β must be dimensionless"))
    argcheck_tensortrace(C,   A, p, q)
    dimcheck_tensortrace(C, α*A, p, q)
    return C
end

function tensorcontract!(C::AbstractAxesDimensions,
                         A, pA::Index2Tuple, conjA::Bool,
                         B, pB::Index2Tuple, conjB::Bool,
                         pAB::Index2Tuple,
                         α::AbstractDimensions, β::AbstractDimensions,
                         backend, allocator)
    β == one(β) || throw(DimensionMismatch("β must be dimensionless"))
    argcheck_tensorcontract(C,   A, pA, B, pB, pAB)
    dimcheck_tensorcontract(C, α*A, pA, B, pB, pAB)
    return C
end

############################# Dimensional checks ##############################

function dimcheck_tensoradd(C::AbstractAxesDimensions,
                            A::AbstractArray, pA::Index2Tuple)
    dimsA, dimsC = normdims.((A, C))
    scaleA, scaleC = dimscale.((A, C))
    scaleC == scaleA ||
        throw(DimensionMismatch("non-matching dimscales"))
    TupleTools.getindices(dimsA, linearize(pA)) == dimsC ||
        throw(DimensionMismatch("non-matching dimensions in uncontracted indices"))
    return nothing
end

function dimcheck_tensortrace(C::AbstractAxesDimensions,
                              A::AbstractArray, p::Index2Tuple, q::Index2Tuple)
    dimsA, dimsC = normdims.((A, C))
    scaleA, scaleC = dimscale.((A, C))
    scaleC == scaleA ||
        throw(DimensionMismatch("non-matching dimscales"))
    allmap(match, TupleTools.getindices(dimsA, q[1]), TupleTools.getindices(dimsA, q[2])) ||
        throw(DimensionMismatch("non-matching dimensions in traced indices"))
    TupleTools.getindices(dimsA, linearize(p)) == dimsC ||
        throw(DimensionMismatch("non-matching dimensions in uncontracted indices"))
    return nothing
end

function dimcheck_tensorcontract(C::AbstractAxesDimensions,
                                 A::AbstractArray, pA::Index2Tuple,
                                 B::AbstractArray, pB::Index2Tuple,
                                 pAB::Index2Tuple)
    dimsA, dimsB, dimsC = normdims.((A, B, C))
    scaleA, scaleB, scaleC = dimscale.((A, B, C))
    scaleC == scaleA * scaleB ||
        throw(DimensionMismatch("non-matching dimscales"))
    allmap(match, TupleTools.getindices(dimsA, pA[2]), TupleTools.getindices(dimsB, pB[1])) ||
        throw(DimensionMismatch("non-matching dimensions in contracted indices"))
    dimsAB = (TupleTools.getindices(dimsA, pA[1])..., TupleTools.getindices(dimsB, pB[2])...)
    TupleTools.getindices(dimsAB, linearize(pAB)) == dimsC ||
        throw(DimensionMismatch("non-matching dimensions in uncontracted indices"))
    return nothing
end

end # of module