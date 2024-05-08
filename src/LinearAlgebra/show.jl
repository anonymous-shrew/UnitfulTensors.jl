function replace_in_print_matrix(A::AbstractUnitfulTensor, i::Integer, j::Integer, s::AbstractString)
    return replace_in_print_matrix(values(A), i, j, s)
end

function show(io::IO, ::MIME{Symbol("text/plain")}, Q::UnitfulTensor{<:Any, <:Any, <:Any, <:AbstractQ})
    print(io, summary(Q))
end

function show(io::IO, mime::MIME{Symbol("text/plain")}, F::AbstractUnitfulFactorization)
    if _issuccess(F)
        summary(io, F)
        for property in _displayed_properties(F)
            _show_property(io, mime, F, property)
        end
    else
        print(io, "Failed factorization of type $(typeof(F))")
    end
end

_issuccess(F::AbstractUnitfulFactorization) = _issuccess(values(F))
_issuccess(F::Factorization) = true

function _show_property(io::IO, mime::MIME{Symbol("text/plain")}, F::AbstractUnitfulFactorization, p::Symbol)
    _show_property(io, mime, F, Val(p))
end
    
function _show_property(io::IO, mime::MIME{Symbol("text/plain")}, F::AbstractUnitfulFactorization, ::Val{p}) where p
    println(io, "\n", _property_description(F, p), ":")
    show(io, mime, getproperty(F, p))
end

function _property_description(F::AbstractUnitfulFactorization, p::Symbol)
    descriptions = _property_descriptions(F)
    if p ∈ keys(descriptions)
        return descriptions[p]
    elseif p === :p
        return "permutation"
    else
        return string(p, " factor")
    end
end

_property_descriptions(F::AbstractUnitfulFactorization) = _property_descriptions(values(F))
_property_descriptions(F::Factorization) = (; )