using ADTI
using CUDA
using LinearAlgebra
using Random
using TeneT
using Zygote
using TeneT:_fit_spaces
using TensorKit
using OptimKit
using ADTI: build_M
#####################################    parameters      ###################################
Random.seed!(42)
atype = Array
pattern = [1;;]
d = ℤ₂Space(0=>2, 1=>2)
D = ℤ₂Space(0=>1, 1=>1)
χ = ℤ₂Space(0=>10, 1=>10)
No = 0
λ = 0
model = Hubbard(1.0, 0.0, 0.0)
system = :fermion
lattice = :square
folder = "data/$model/$system/$lattice/$pattern/"
boundary_alg = VUMPS(ifupdown=true,
                     ifdownfromup=true, 
                     ifsimple_eig=true,
                     maxiter=30, 
                     miniter=1,
                     maxiter_ad=10,
                     miniter_ad=3,  
                     ifcheckpoint=false,
                     verbosity=2,
                     show_every=10
)
params = iPEPSOptimize{system, lattice}(boundary_alg=boundary_alg,
                                        optimizer=LBFGS(200; gradtol=1e-10, maxiter=100),
                                        reuse_env=true, 
                                        ifcheckpoint=false, 
                                        ifflatten=true,
                                        ifNN=false,
                                        verbosity=4, 
                                        λ=λ,
                                        folder=folder
)
A = init_ipeps(;atype, params, No, d, D, χ, pattern)
# M = build_M(A, params)
# @show M
# A = ADTI.init_ipeps_spin111(;atype, model, params, No, ifWp=true, ϵ = 0, χ, Ni, Nj)
# A = ADTI.init_ipeps_h5(;atype, params, ifWp=false, ifreal=true, model, file="./data/kitsShf_sikh3nfr1D8D8.h5", D, Ni, Nj)
############################################################################################

function _restriction_ipeps(A)
    # Ar = Zygote.Buffer(A)
    # len = length(unique(pattern))
    # for i in 1:len
    #     temp = A[i] + _fit_spaces(permute(A[i], ((2,1,3,4,5),())), A[i])
    #     temp = temp + _fit_spaces(permute(temp, ((1,2,3,5,4),())), temp)
    #     # temp = temp + _fit_spaces(permute(temp, ((1,5,3,4,2),())), temp)
    #     # temp = temp + _fit_spaces(permute(temp, ((2,1,3,5,4),())), temp)
    #     Ar[i] = temp
    # end
    # Ar = copy(Ar)
    # return Ar/norm(Ar)
    return A
end
optimise_ipeps(A, model, χ, params;
               restriction_ipeps = _restriction_ipeps);