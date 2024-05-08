# API

------------

# Scalar API

## ​ Types

```@docs
AbstractUnitfulScalar
UnitfulScalar
AbstractDimensions
SIDimensions
```

## ​ Constants

```@docs
NoDims
```

## ​ Functions

```@docs
value
values(::AbstractUnitfulScalar)
dimensions(::AbstractUnitfulScalar)
dimexps
```

# Tensor API
## ​ Basic types

```@docs
AbstractUnitfulTensor
UnitfulTensor
AbstractAxesDimensions
AxesDimensions
AbstractAxisDimensions
AxisDimensions
```

## ​ Extra types

```@docs
AdjointAxesDimensions
AbstractUnitfulFactorization
UnitfulFactorization
UnitfulGeneralizedFactorization
```

## ​ Type aliases

```@docs
AbstractUnitfulVector
UnitfulVector
AbstractUnitfulMatrix
UnitfulMatrix
AbstractVectorDimensions
AbstractMatrixDimensions
AdjointVectorDimensions
AdjointMatrixDimensions
```

## ​ Constructors

```@docs
UnitfulTensor(::AbstractArray)
AxesDimensions(A::AbstractArray{<:AbstractDimensions})
AxesDimensions(dims, scale; normalize = true)
AxisDimensions(dims::AbstractVector{<:AbstractDimensions}; normalize = true)
```

## ​ Functions

```@docs
values(::AbstractUnitfulTensor)
dimensions(::AbstractUnitfulTensor)
normdims
dimscale
dimsplat
dimsvec
nodims
tensor_product
tensor_factorize
ishomogeneous
issquareable
isendomorphic
match(::AbstractAxisDimensions, ::AbstractAxisDimensions)
units_off
tensorize_literals
```
