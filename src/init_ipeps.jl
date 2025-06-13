"""
    init_ipeps(;atype = Array, Ni::Int, Nj::Int, D::Int)

return a random `ipeps` with bond dimension `D` and physical dimension 2.
"""
function init_ipeps(;atype = Array, params::iPEPSOptimize{S,L}, No::Int=0, d::Int, D::Int, χ::Int) where {S,L}
    if No != 0
        file = joinpath(params.folder, "D$(D)_χ$(χ)/ipeps/ipeps_No.$(No).jld2")
        @info "load ipeps from file: $file"
        A = load(file, "bcipeps")
    else
        N = length(unique(params.pattern))
        Ni, Nj = size(params.pattern)
        @info "generate random $S $L ipeps"
        if S == :boson
            if L == :merge
                A = rand(ComplexF64, D,D,D,D,d^2, N) .- 0.5
            elseif L == :brickwall
                Ni % 2 == 0 && Nj % 2 == 0 || throw(ArgumentError("Ni and Nj should be even"))
                A = rand(ComplexF64, D,1,D,D,d, N) .- 0.5
            end
        elseif S == :fermion
            A = rand(ComplexF64, D,D,d,D,D, N) .- 0.5
        end
        A /= norm(A)
    end
    return atype(A)
end

function init_ipeps_spin111(;atype = Array, model, params::iPEPSOptimize{F}, No::Int=0, ifWp::Bool=false, ϵ::Real=1e-3, χ::Int, Ni::Int, Nj::Int) where F
    d = Int(2*model.S + 1) 
    if No != 0
        file = joinpath(params.folder, "D$(D)_χ$(χ)/ipeps/ipeps_No.$(No).jld2")
        @info "load ipeps from file: $file"
        A = load(file, "bcipeps")
    else
        @info "generate spin-(1,1,1) ipeps"
        spin111 = [1,sqrt(2)/2*(1-1im),1]
        spin111 = spin111 / norm(spin111)
        if F == :merge
            A = reshape(ein"a,b->ab"(spin111, spin111), 1,1,1,1,d^2, Ni,Nj)
        elseif F == :brickwall
            A = rand(ComplexF64, D,1,D,D,d, Ni,Nj) .- 0.5
        end
        A /= norm(A)
    end
    if ifWp
        Wp = _arraytype(A)(build_Wp(model.S, params))
        A = build_A(A + randn(ComplexF64, size(A)) * ϵ, Wp, params) 
    end
    return atype(A)
end

function init_ipeps_h5(;atype = Array, params, model, ϵ=0, ifWp=false, ifreal=false, file, D::Int, Ni::Int, Nj::Int)
    d = Int(2*model.S + 1) 
    A = zeros(ComplexF64, D,1,D,D,d, Ni,Nj)
    @info "load SU init ipeps from $file"
    for i in 1:Ni, j in 1:Nj
        if (i+j) % 2 == 0
            if ifreal
                A[:,:,:,:,:,i,j] = permutedims(h5open(file, "r") do file
                    read(file, "$i$(7-j)r")
                end
                , (5,2,3,4,1))
            else
                A[:,:,:,:,:,i,j] = permutedims(h5open(file, "r") do file
                    read(file, "$i$(7-j)r")
                end + 1im * h5open(file, "r") do file
                    read(file, "$i$(7-j)i")
                end
                , (5,2,3,4,1))
            end
        else
            if ifreal
                A[:,:,:,:,:,i,j] = permutedims(h5open(file, "r") do file
                    read(file, "$i$(7-j)r")
                end, (3,4,5,2,1))
            else
                A[:,:,:,:,:,i,j] = permutedims(h5open(file, "r") do file
                    read(file, "$i$(7-j)r")
                end + 1im * h5open(file, "r") do file
                    read(file, "$i$(7-j)i")
                end, (3,4,5,2,1))
            end
        end
    end
    if ifWp
        Wp = _arraytype(A)(build_Wp(model.S, params))
        A = build_A(A + randn(ComplexF64, size(A)) * ϵ, Wp, params) 
    end
    return atype(A)
end

function _restriction_ipeps(A)
    return A/norm(A)
end

function init_ipeps_form_small_spin(;atype = Array, file, D::Int, d::Int, Ni::Int, Nj::Int)
    @info "load ipeps from file: $file"
    Ao = load(file, "bcipeps")
    A = zeros(ComplexF64, D,D,D,D,d,Ni,Nj)
    A[:,:,:,:,[6,7,10,11],:,:] = Ao
    return atype(A)
end