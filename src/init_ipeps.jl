"""
    init_ipeps(;atype = Array, Ni::Int, Nj::Int, D::Int)

return a random `ipeps` with bond dimension `D` and physical dimension 2.
"""
function init_ipeps(;atype = Array, No::Int=0, D, χ, Ni::Int, Nj::Int, params::iPEPSOptimize{:boson,L}) where {L}
    if No != 0
        file = joinpath(params.folder, "D$(D)_χ$(χ)/ipeps/ipeps_No.$(No).jld2")
        @info "load ipeps from file: $file"
        A = load(file, "bcipeps")
    else
        @info "generate random boson $L ipeps"
        if L == :merge
            A = rand(ComplexF64, D,D,D,D,d^2, Ni,Nj) .- 0.5
        elseif L == :brickwall
            Ni % 2 == 0 && Nj % 2 == 0 || throw(ArgumentError("Ni and Nj should be even"))
            A = rand(ComplexF64, D,1,D,D,d, Ni,Nj) .- 0.5
        end
        A /= norm(A)
    end
    return atype(A)
end

function init_ipeps(;atype = Array, No::Int=0, d, Ni::Int, Nj::Int, params::iPEPSOptimize{:fermion,L}) where {L}
    if No != 0
        file = joinpath(params.folder, "D$(D)_χ$(χ)/ipeps/ipeps_No.$(No).jld2")
        @info "load ipeps from file: $file"
        A = load(file, "bcipeps")
    else
        @info "generate random fermion $L ipeps"
        sitetype = params.sitetype
        if sitetype === nothing
            A = rand(ComplexF64, D,D,d,D,D, Ni,Nj) .- 0.5
        else
            qnD, _, dimsD, _ = params.boundary_alg.U1info
            D = sum(dimsD)
            randdims = sum(prod.(
					            zerosinitial(Val(:U1), atype, ComplexF64, D, D, d, D, D; 
								dir = [-1, -1, 1, 1, 1], 
								indqn = [qnD, qnD, getqrange(sitetype, d)..., qnD, 	qnD],                    
								indims = [dimsD, dimsD, getblockdims(sitetype, d)..., dimsD, dimsD], 
								f=[0],
								ifZ2=sitetype.ifZ2
								).dims))
            A = rand(ComplexF64, randdims, Ni,Nj) .- 0.5
        end
        A /= norm(A)
    end
    return atype(A)
end

function _restriction_ipeps(A)
    return A/norm(A)
end