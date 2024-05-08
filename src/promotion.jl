# A promotion system that enables mixing UnitfulTensors and Arrays of Numbers in arithmetic operations

promote_rule(::Type{T}, ::Type{<:Number}) where T <: AbstractUnitfulScalar = T
(::Type{T})(x::Number) where T <: UnitfulScalar = T(x, NoDims)

promote_unitful(f::Union{Function, Type},
    vals, dims::AbstractDimensions, args...) = _promote_unitful_scalar(args...)
promote_unitful(f::Union{Function, Type},
    vals, dims::AbstractAxesDimensions, args...) = _promote_unitful_tensor(args...)
promote_unitful(f::Union{Function, Type},
    vals, dims::Tuple{Vararg{AbstractAxesDimensions}}, args...) = _promote_unitful_tensor(args...)

_scalar_type(x::T) where T = _scalar_type(T)
_scalar_type(x::AbstractArray) = _scalar_type(eltype(x))
_scalar_type(x::Type) = throw(ArgumentError("Unsupported type $x in _scalar_type"))
_scalar_type(::Type{T}) where T <: Number = T
_scalar_type(::Type{<:UnitfulScalar}) = UnitfulScalar

_promote_unitful_scalar(xs...) = _promote_unitful_scalar(_scalar_type.(xs)...)
_promote_unitful_scalar(xs::Type{<:Number}...) = afoldl(_promote_unitful_scalar, xs...) # foldl is slower
_promote_unitful_scalar() = UnitfulScalar
_promote_unitful_scalar(::Type{<:Number}) = UnitfulScalar
_promote_unitful_scalar(::Type{T}) where T <: AbstractUnitfulScalar = T
_promote_unitful_scalar(::Type{<:Number}, ::Type{<:Number}) = UnitfulScalar
_promote_unitful_scalar(::Type{<:AbstractUnitfulScalar}, ::Type{<:AbstractUnitfulScalar}) = UnitfulScalar
_promote_unitful_scalar(::Type{T}, ::Type{T}) where T <: AbstractUnitfulScalar = T
_promote_unitful_scalar(::Type{T}, ::Type{<:Number}) where T <: AbstractUnitfulScalar = T
_promote_unitful_scalar(::Type{<:Number}, ::Type{T}) where T <: AbstractUnitfulScalar = T

_tensor_type(x) = UnitfulTensor
_tensor_type(x::UnitfulTensor) = UnitfulTensor
_tensor_type(x::T) where T <: AbstractArray = T

_promote_unitful_tensor(xs...) = _promote_unitful_tensor(_tensor_type.(xs)...)
_promote_unitful_tensor(xs::Type{<:AbstractArray}...) = afoldl(_promote_unitful_tensor, xs...) # foldl is slower
_promote_unitful_tensor() = UnitfulTensor
_promote_unitful_tensor(::Type{<:AbstractArray}) = UnitfulTensor
_promote_unitful_tensor(::Type{T}) where T <: AbstractUnitfulTensor = T
_promote_unitful_tensor(::Type{<:AbstractArray}, ::Type{<:AbstractArray}) = UnitfulTensor
_promote_unitful_tensor(::Type{<:AbstractUnitfulTensor}, ::Type{<:AbstractUnitfulTensor}) = UnitfulTensor
_promote_unitful_tensor(::Type{T}, ::Type{T}) where T <: AbstractUnitfulTensor = T
_promote_unitful_tensor(::Type{T}, ::Type{<:AbstractArray}) where T <: AbstractUnitfulTensor = T
_promote_unitful_tensor(::Type{<:AbstractArray}, ::Type{T}) where T <: AbstractUnitfulTensor = T

promote_dims(f::Union{Function, Type}, dims, args...) = _promote_axesdims(args...)

_axes_type(x) = AxesDimensions

_promote_axesdims(xs...) = _promote_axesdims(_axes_type.(xs)...)
_promote_axesdims(xs::Type{<:AbstractAxesDimensions}...) = afoldl(_promote_axesdims, xs...) # foldl is slower
_promote_axesdims() = AxesDimensions
_promote_axesdims(::Type{T}) where T <: AbstractAxesDimensions = T
_promote_axesdims(::Type{<:AbstractAxesDimensions}, ::Type{<:AbstractAxesDimensions}) = AxesDimensions
_promote_axesdims(::Type{T}, ::Type{T}) where T <: AbstractAxesDimensions = T