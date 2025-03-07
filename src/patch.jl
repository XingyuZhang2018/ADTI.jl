rejoinpath(S::String...) = replace(Base.joinpath(S), "=>" => "→")
@non_differentiable rejoinpath(S...)

function Base.:+(a::NamedTuple, b::StructArray)
    return StructArray(a.data + b.data, b.pattern)
end