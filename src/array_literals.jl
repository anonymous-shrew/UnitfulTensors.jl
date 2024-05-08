# Not sure if this is a sane idea
const _tensorize_literals_body = quote
    struct Arrayable{T} x::T end
    
    for f in (:value, :values, :dimensions) eval(:($f(x::Arrayable) = $f(x.x))) end
    
    convert(::Type{Arrayable{T}}, x) where T = Arrayable(convert(T, x))
    convert(::Type{T}, x::T) where T <: Arrayable = x
    
    for f in (:vect, :hcat, :vcat)
        eval(:(
                $f(xs::T...) where T <: AbstractUnitfulScalar = UnitfulTensor($f(Arrayable{T}.(xs)...))
                ))
    end
    function hvcat(rows::Tuple{Vararg{Int}}, xs::T...) where T <: AbstractUnitfulScalar
        return UnitfulTensor(hvcat(rows, Arrayable{T}.(xs)...))
    end
    
    getindex(::Type{T}, xs...) where T <: AbstractUnitfulScalar = UnitfulTensor(Arrayable{T}[xs...])
    
    function hvcat_fill!(A::Array{<:AbstractUnitfulScalar}, xs::Tuple)
        B = similar(Array{Arrayable{eltype(A)}}, axes(A))
        return UnitfulTensor(hvcat_fill!(B, xs))
    end
    
    function typed_hvcat(::Type{T}, rows::Tuple{Vararg{Int}}, xs::Number...) where T <: AbstractUnitfulScalar
        return UnitfulTensor(typed_hvcat(Arrayable{T}, rows, Arrayable{T}.(xs)...))
    end
end

"""
    tensorize_literals()

Make array literals involving [`AbstractUnitfulScalar`](@ref)s return [`UnitfulTensor`](@ref)s
instead of `Array`s, so that you don't have to write `UnitfulTensor([...])` all the time.

This feature is experimental. It might break something and may be removed in the future.
"""
tensorize_literals() = eval(_tensorize_literals_body)
