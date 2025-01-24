using ADIT
using CUDA
using LinearAlgebra
using Random
using TeneT
using Zygote

#####################################    parameters      ###################################
Random.seed!(100)
atype = Array
Ni, Nj = 1, 1
d, D, χ = 4, 3, 20
No = 0
model = Hubbard(1.0, 0.0, 0.0)
system = :fermion
lattice = :square
folder = "data/$model/$system/$lattice/$(Ni)x$(Nj)/"
boundary_alg = VUMPS(ifupdown=true,
                     ifdownfromup=false, 
                     ifsimple_eig=true,
                     maxiter=30, 
                     miniter=1,
                     maxiter_ad=3,
                     miniter_ad=3,
                     ifcheckpoint=true,
                     verbosity=3,
                     show_every=10
)
params = iPEPSOptimize{system, lattice}(boundary_alg=boundary_alg,
                                        reuse_env=true, 
                                        ifcheckpoint=true, 
                                        ifflatten=true,
                                        ifNN=false,
                                        verbosity=4, 
                                        maxiter=1000,
                                        tol=1e-10,
                                        folder=folder
)
A = init_ipeps(;atype, model, params, No, d, D, χ, Ni, Nj)
# A = ADIT.init_ipeps_spin111(;atype, model, params, No, ifWp=true, ϵ = 0, χ, Ni, Nj)
# A = ADIT.init_ipeps_h5(;atype, params, ifWp=false, ifreal=true, model, file="./data/kitsShf_sikh3nfr1D8D8.h5", D, Ni, Nj)
############################################################################################

function _restriction_ipeps(A)
    Ar = Zygote.Buffer(A)
    Ni, Nj = size(A)[[6,7]]
    for j in 1:Nj, i in 1:Ni
        Ar[:,:,:,:,:,i,j] = A[:,:,:,:,:,i,j] + permutedims(A[:,:,:,:,:,i,j], (2,1,3,5,4))
    end
    Ar = copy(Ar)
    return Ar/norm(Ar)
    # return A/norm(A)
end
optimise_ipeps(A, model, χ, params;
               restriction_ipeps = _restriction_ipeps)