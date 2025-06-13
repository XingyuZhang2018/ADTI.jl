using ADTI
using CUDA
using LinearAlgebra
using Random
using TeneT
using Zygote
using OptimKit

#####################################    parameters      ###################################
Random.seed!(100)
atype = Array
pattern = [1 2; 2 1]
Ni, Nj = size(pattern)
d, D, χ = 4, 2, 20
No = 10
model = Hubbard(1.0, 0.0, 0.0)
system = :fermion
lattice = :square
folder = "data/$model/$system/$lattice/$(pattern)/"
boundary_alg = VUMPS(ifupdown=true,
                     ifdownfromup=false, 
                     ifparallelupdown=false,
                     ifsimple_eig=true,
                     maxiter=300, 
                     miniter=1,
                     maxiter_ad=10,
                     miniter_ad=3,
                     ifcheckpoint=false,
                     power_iter=5,
                     power_iter_obs=50,
                     verbosity=3,
                     show_every=10
)
params = iPEPSOptimize{system, lattice}(pattern=pattern,
                                        boundary_alg=boundary_alg,
                                        optimizer=LBFGS(200; maxiter=100, verbosity=1, gradtol=1e-7),
                                        reuse_env=true, 
                                        ifcheckpoint=false, 
                                        ifflatten=true,
                                        ifNN=false,
                                        verbosity=4, 
                                        maxiter=0,
                                        tol=1e-10,
                                        folder=folder
)
A = init_ipeps(;atype, params, No, d, D, χ)
# A = ADTI.init_ipeps_spin111(;atype, model, params, No, ifWp=true, ϵ = 0, χ, Ni, Nj)
# A = ADTI.init_ipeps_h5(;atype, params, ifWp=false, ifreal=true, model, file="./data/kitsShf_sikh3nfr1D8D8.h5", D, Ni, Nj)
############################################################################################

function _restriction_ipeps(A)
    return A
end
optimise_ipeps(A, model, χ, params;
               restriction_ipeps = _restriction_ipeps)