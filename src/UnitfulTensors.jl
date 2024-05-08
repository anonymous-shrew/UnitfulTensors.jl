module UnitfulTensors

include("FastQuantities.jl")

using .FastQuantities
using .FastQuantities: allmap
import .FastQuantities: value, dimensions # types.jl, array_literals.jl

using Base: delete_method, # types.jl
            afoldl, # promotion.jl
            tail, # tensor_product.jl
            index_ndims, index_shape, _unsafe_getindex!, __maybe_reshape, LogicalIndex, Slice # indexing.jl

using Base.IteratorsMD: flatten, split # tensor_product.jl

import Base: values, match, ==, # types.jl
             promote_rule, # promotion.jl
             IndexStyle, size, getindex, setindex!, _unsafe_getindex, reshape, _maybe_reshape, first, # indexing.jl
             parent, # adjtrans.jl
             ≈, +, -, *, \, /, ^, sqrt, # arithmetic.jl           
             exp ,  log ,  cis , # transcendental.jl
             sin ,  cos ,  tan ,  cot ,  sec ,  csc ,
            asin , acos , atan , acot , asec , acsc ,
             sinh,  cosh,  tanh,  coth,  sech,  csch,
            asinh, acosh, atanh, acoth, asech, acsch,
             sind,  cosd,  tand,# coth,  sech,  csch, (they don't yet accept matrices in stdlib)
            asind, acosd, atand, acotd, asecd, acscd,
            sincos, sincosd,
             iterate, getproperty, # factorization.jl
             show, replace_in_print_matrix, # show.jl
             convert, vect, hcat, vcat, hvcat, hvcat_fill!, typed_hvcat # array_literals.jl

using LinearAlgebra: AbstractQ, # types.jl
                     Adjoint, Transpose, # adjtrans.jl
                     AdjointAbsVec, TransposeAbsVec, AdjOrTransAbsVec, ⋅, # arithmetic.jl
                     Factorization, sym_uplo, # factorization.jl 
                     LU, # factorizations.jl
                     BunchKaufman, Cholesky, CholeskyPivoted, LDLt,
                     Eigen, Hessenberg, Schur, SVD,
                     GeneralizedEigen, GeneralizedSchur, GeneralizedSVD,
                     LQ, QR, QRCompactWY, QRPivoted,
                     NoPivot, RowMaximum, RowNonZero, ColumnNorm,
                     checksquare

import LinearAlgebra: adjoint, adjoint!, transpose, transpose!, # adjtrans.jl
                      issymmetric, ishermitian, isposdef, isposdef!, # is.jl
                      isdiag, istril, istriu,
                      dot, cross, kron, kron!, lmul!, rmul!, ldiv!, rdiv!, # arithmetic.jl
                      inv, pinv,
                      det, logdet, logabsdet, tr, # generic.jl
                      lyap, sylvester,
                      rotate!, reflect!, givens,
                      cond, condskeel, norm, normalize, normalize!, nullspace, opnorm, rank, # norm.jl
                      issuccess, # factorization.jl
                      lu, lu!, # factorizations.jl
                      bunchkaufman, bunchkaufman!, cholesky, cholesky!, ldlt, ldlt!,
                      eigen, eigen!, hessenberg, hessenberg!, schur, schur!, svd, svd!,
                      lq, lq!, qr, qr!,
                      eigvals, eigvals!, eigvecs, eigmax, eigmin, svdvals, svdvals!,
                      lowrankupdate, lowrankupdate!, lowrankdowndate, lowrankdowndate!,
                      diag, hermitianpart, hermitianpart!, tril, tril!, triu, triu!, # structured.jl
                      Diagonal

import SparseArrays: sparse # sparse.jl
       
export # FastQuantities.jl
       AbstractDimensions, SIDimensions, 
       AbstractUnitfulScalar, UnitfulScalar,
       NoDims, 𝐓, 𝐋, 𝐌, 𝐈, 𝚯, 𝐍, 𝐉,
       dimexps, value, values, dimensions,
       @u_str,

       # types.jl
       AbstractAxisDimensions, AxisDimensions,
       AbstractAxesDimensions, AxesDimensions,
       AbstractVectorDimensions, AbstractMatrixDimensions,
       AbstractVecOrMatDimensions, AbstractDimensionsOrAxesDimensions,
       AbstractUnitfulTensor, UnitfulTensor,
       AbstractUnitfulVector, UnitfulVector,
       AbstractUnitfulMatrix, UnitfulMatrix,
       AbstractUnitfulVecOrMat, UnitfulVecOrMat,
       AbstractUnitfulScalarOrTensor, UnitfulScalarOrTensor,

       # adjtrans.jl
       AdjointAxesDimensions, AdjointVectorDimensions, AdjointMatrixDimensions,

       # factorization.jl
       AbstractUnitfulFactorization, UnitfulFactorization, UnitfulGeneralizedFactorization,

       # types.jl
       dimsvec, normdims, dimscale, dimsplat,
       nodims,

       # tensor_product.jl
       tensor_factorize, tensor_product, ⊗,

       # is.jl
       ishomogeneous, issquareable, isendomorphic,

       # types.jl, array_literals.jl
       units_off, tensorize_literals

include("types.jl")
include("promotion.jl")
include("tensor_product.jl")
include("indexing.jl")

include("LinearAlgebra/adjtrans.jl")
include("LinearAlgebra/is.jl")
include("LinearAlgebra/arithmetic.jl")
include("LinearAlgebra/transcendental.jl")
include("LinearAlgebra/generic.jl")
include("LinearAlgebra/norm.jl")

include("LinearAlgebra/factorization.jl")
include("LinearAlgebra/show.jl")
include("LinearAlgebra/factorizations.jl")

include("LinearAlgebra/structured.jl")

include("LinearAlgebra/sparse.jl")

include("array_literals.jl")

end # of module
