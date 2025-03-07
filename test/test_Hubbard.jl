using ADTI
using CUDA
using LinearAlgebra
using Random
using TeneT
using Zygote
using TeneT:_fit_spaces
using TensorKit
using OptimKit

#####################################    parameters      ###################################
Random.seed!(100)
atype = Array
pattern = [1 2; 2 1]
d = ℤ₂Space(0=>2, 1=>2)
D = ℤ₂Space(0=>1, 1=>1)
χ = ℤ₂Space(0=>10, 1=>10)
No = 0
model = Hubbard(1.0, 12.0, 6.0)
system = :fermion
lattice = :square
folder = "data/$model/$system/$lattice/$pattern/"
boundary_alg = VUMPS(ifupdown=true,
                     ifdownfromup=false, 
                     ifsimple_eig=true,
                     maxiter=30, 
                     miniter=1,
                     maxiter_ad=3,
                     miniter_ad=3,
                     ifcheckpoint=false,
                     verbosity=3,
                     show_every=10
)
params = iPEPSOptimize{system, lattice}(boundary_alg=boundary_alg,
                                        optimizer=LBFGS(; gradtol=1e-10, maxiter=100),
                                        reuse_env=true, 
                                        ifcheckpoint=false, 
                                        ifflatten=true,
                                        ifNN=false,
                                        verbosity=4, 
                                        folder=folder
)
A = init_ipeps(;atype, params, No, d, D, χ, pattern)
# A = ADTI.init_ipeps_spin111(;atype, model, params, No, ifWp=true, ϵ = 0, χ, Ni, Nj)
# A = ADTI.init_ipeps_h5(;atype, params, ifWp=false, ifreal=true, model, file="./data/kitsShf_sikh3nfr1D8D8.h5", D, Ni, Nj)
############################################################################################

function _restriction_ipeps(A)
    # Ar = Zygote.Buffer(A)
    # len = length(unique(pattern))
    # for i in 1:len
    #     Ar[i] = A[i] + _fit_spaces(permute(A[i], ((2,1,3,5,4),())), A[i])
    # end
    # Ar = copy(Ar)
    # return Ar/norm(Ar)
    return A/norm(A)
end
optimise_ipeps(A, model, χ, params;
               restriction_ipeps = _restriction_ipeps);