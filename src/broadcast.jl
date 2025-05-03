struct UnitfulTensorStyle{N} <: AbstractArrayStyle{N} end
struct AxesDimensionsStyle{N} <: AbstractArrayStyle{N} end

UnitfulTensorStyle{M}(::Val{N}) where {N,M} = UnitfulTensorStyle{N}()
AxesDimensionsStyle{M}(::Val{N}) where {N,M} = AxesDimensionsStyle{N}()

BroadcastStyle(T::Type{<:AbstractUnitfulTensor}) = UnitfulTensorStyle{ndims(T)}()
BroadcastStyle(T::Type{<:AbstractUnitfulScalar}) = UnitfulTensorStyle{0}()
BroadcastStyle(T::Type{<:AbstractAxesDimensions}) = AxesDimensionsStyle{ndims(T)}()
BroadcastStyle(T::Type{<:AbstractDimensions}) = AxesDimensionsStyle{0}()

similar(b::Broadcasted{<:AxesDimensionsStyle}, shape) = similar(Array{combine_eltypes(b.f, b.args)}, shape)

function copy(x::Broadcasted{UnitfulTensorStyle{N}}) where N
    dims = copy(convert(Broadcasted{AxesDimensionsStyle{N}}, broadcasted(dimensions, x)))
    vals = broadcast(values ∘ x.f, UnitEater.(x.args)...)
    T = (N == 0) ? UnitfulScalar : UnitfulTensor
    return T(vals, dims)
end

copy(x::Broadcasted{<:AxesDimensionsStyle}) = _get_single_index(x, eachindex(x); unsafe=true)
copy(x::Broadcasted{AxesDimensionsStyle{0}}) = x[]

function copyto!(dest::AbstractUnitfulTensor, x::Broadcasted{UnitfulTensorStyle{N}}) where N
    copyto!(dimensions(dest), convert(Broadcasted{AxesDimensionsStyle{N}}, broadcasted(dimensions, x)))
    broadcast!(values ∘ x.f, values(dest), UnitEater.(x.args)...)
    return dest
end

function copyto!(dest::AbstractAxesDimensions, x::Broadcasted{<:AxesDimensionsStyle})
    result = _get_single_index(x, eachindex(x); unsafe=true)
    if dest != result
        throw(DimensionMismatch("destination has dimensions $dest; result has dimensions $result"))
    end
    return dest
end


UnitEater(x::AbstractUnitfulScalar) = UnitEater(values(x))
UnitEater(x::AbstractUnitfulTensor) = UnitEater(values(x))
UnitEater(x::Number) = ScalarUnitEater(x)
UnitEater(x::AbstractArray) = TensorUnitEater(x)
UnitEater(x::Broadcasted) = broadcasted(x.f, UnitEater.(x.args)...)
UnitEater(x) = x

struct ScalarUnitEater{T} <: Number x::T end
struct TensorUnitEater{T, N, AT <: AbstractArray{T, N}} <: AbstractArray{ScalarUnitEater{T}, N} x::AT end

values(x::ScalarUnitEater) = x.x

size(x::TensorUnitEater) = size(x.x)
axes(x::TensorUnitEater) = axes(x.x)
IndexStyle(x::TensorUnitEater) = IndexStyle(x.x)
getindex(x::TensorUnitEater, I::Vararg{Int, N}) where N = UnitEater(x.x[I...])

promote_rule(::Type{T}, ::Type) where T <: ScalarUnitEater = T
promote_rule(::Type{ScalarUnitEater{T1}}, ::Type{ScalarUnitEater{T2}}) where {T1, T2} = ScalarUnitEater{promote_type(T1, T2)}
convert(::Type{T}, x::Number) where T <: ScalarUnitEater = T(values(x))
convert(::Type{T}, x::T) where T <: ScalarUnitEater = x

_promote_unitful_scalar(::Type{<:ScalarUnitEater{T1}}, ::Type{<:AbstractUnitfulScalar{T2}}) where {T1, T2} = ScalarUnitEater{promote_type(T1, T2)}
_promote_unitful_scalar(::Type{<:AbstractUnitfulScalar{T1}}, ::Type{<:ScalarUnitEater{T2}}) where {T1, T2} = ScalarUnitEater{promote_type(T1, T2)}

for f in (:+, :-, :*, :/, :^)
    # TODO: replace with promotion
    @eval $f(x::ScalarUnitEater, y::AbstractUnitfulScalar) = ScalarUnitEater($f(values.((x, y))...))
    @eval $f(x::AbstractUnitfulScalar, y::ScalarUnitEater) = ScalarUnitEater($f(values.((x, y))...))
    @eval $f(x::T, y::T) where {T<:ScalarUnitEater} = ScalarUnitEater($f(values.((x, y))...))
end