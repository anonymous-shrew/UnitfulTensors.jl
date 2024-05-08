# https://github.com/JuliaLang/julia/issues/53380
@test SIDimensions <: AbstractDimensions
@test AxisDimensions <: AbstractAxisDimensions
@test AxesDimensions <: AbstractAxesDimensions
@test UnitfulScalar <: AbstractUnitfulScalar
@test UnitfulTensor <: AbstractUnitfulTensor

@test AbstractDimensions <: Number
@test AbstractUnitfulScalar <: Number
@test AbstractAxisDimensions <: AbstractVector
@test AbstractAxesDimensions <: AbstractArray
@test AbstractUnitfulTensor <: AbstractArray

@test AbstractVectorDimensions <: AbstractVecOrMatDimensions
@test AbstractMatrixDimensions <: AbstractVecOrMatDimensions
@test AbstractVecOrMatDimensions <: AbstractAxesDimensions
@test AbstractDimensions <: AbstractDimensionsOrAxesDimensions
@test AbstractAxesDimensions <: AbstractDimensionsOrAxesDimensions

@test AbstractUnitfulVector <: AbstractUnitfulVecOrMat
@test AbstractUnitfulMatrix <: AbstractUnitfulVecOrMat
@test AbstractUnitfulVecOrMat <: AbstractUnitfulTensor
@test AbstractUnitfulScalar <: AbstractUnitfulScalarOrTensor
@test AbstractUnitfulTensor <: AbstractUnitfulScalarOrTensor

@test UnitfulVector <: UnitfulVecOrMat
@test UnitfulMatrix <: UnitfulVecOrMat
@test UnitfulVecOrMat <: UnitfulTensor
@test UnitfulScalar <: UnitfulScalarOrTensor
@test UnitfulTensor <: UnitfulScalarOrTensor

@test AdjointAxesDimensions <: AbstractMatrixDimensions
@test AdjointVectorDimensions <: AdjointAxesDimensions
@test AdjointMatrixDimensions <: AdjointAxesDimensions

@test UnitfulFactorization <: AbstractUnitfulFactorization
@test UnitfulGeneralizedFactorization <: AbstractUnitfulFactorization