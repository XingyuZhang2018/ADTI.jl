"""
    init_ipeps(;atype = Array, Ni::Int, Nj::Int, D::Int)

return a random `ipeps` with bond dimension `D` and physical dimension 2.
"""
function init_ipeps(;atype = Array, params::iPEPSOptimize{S,L}, No::Int=0, d::VectorSpace, D::VectorSpace, χ::VectorSpace, pattern::Matrix{Int}) where {S,L}
    if No != 0
        file = rejoinpath(params.folder, "D$(D)_χ$(χ)/ipeps/ipeps_No.$(No).jld2")
        params.verbosity >= 2 && @info "load ipeps from file: $file"
        A = load(file, "bcipeps")
    else
        params.verbosity >= 2 && @info "generate random $S $L ipeps"
        len = length(unique(pattern))
        Ni, Nj = size(pattern)
        if S == :boson
            if L == :merge
                A = rand(ComplexF64, [D*D*D'*D' ← d^2 for i in 1:len], pattern)
            elseif L == :brickwall
                Ni % 2 == 0 && Nj % 2 == 0 || throw(ArgumentError("Ni and Nj should be even"))
                A = rand(ComplexF64, [D*one(D)*D'*D' ← d for i in 1:len], pattern)
            end
        elseif S == :fermion
            A = rand(ComplexF64, [D'*D'*d*D*D for i in 1:len], pattern)
        end
        A /= norm(A)
    end
    return atype(A)
end

function _restriction_ipeps(A)
    return A/norm(A)
end