using ADTI
using ADTI: swapgate
using CUDA
using LinearAlgebra
using OMEinsum
using Random
using TeneT
using Zygote


#####################################    parameters      ###################################
Random.seed!(100)
atype = Array
Ni, Nj = 1, 1
d, D, χ = 4, 2, 20
No = 1
model = Free_Fermion()
system = :fermion
lattice = :square
folder = "data/$model/$system/$lattice/$(Ni)x$(Nj)/"
boundary_alg = VUMPS(ifupdown=true,
                     ifdownfromup=false, 
                     ifsimple_eig=true,
                     maxiter=100, 
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
# A = ADTI.init_ipeps_spin111(;atype, model, params, No, ifWp=true, ϵ = 0, χ, Ni, Nj)
# A = ADTI.init_ipeps_h5(;atype, params, ifWp=false, ifreal=true, model, file="./data/kitsShf_sikh3nfr1D8D8.h5", D, Ni, Nj)
############################################################################################

function _restriction_ipeps(A)
    Ar = Zygote.Buffer(A)
    Ni, Nj = size(A)[[6,7]]
    SdD = swapgate(d, D)
    restriction1(x) = (x + ein"abcde,cfgb->dfgae"(conj(x), SdD))/2
    restriction2(x) = (x + ein"abcde,cdgf->adgbe"(conj(x), SdD))/2
    restriction3(x) = (x + permutedims(conj(x), (2,1,3,5,4)))/2
    for j in 1:Nj, i in 1:Ni
        # A1 = restriction3(A[:,:,:,:,:,i,j])
        # A2 = restriction1(A1)
        # Ar[:,:,:,:,:,i,j] = A2
        Ar[:,:,:,:,:,i,j] = restriction3(A[:,:,:,:,:,i,j])
        # @show norm(Ar[:,:,:,:,:,i,j] - ein"abcde,cfgb->dfgae"(conj(Ar[:,:,:,:,:,i,j]), SdD))
    end
    Ar = copy(Ar)
    return Ar/norm(Ar)
    # return A/norm(A)
end
optimise_ipeps(A, model, χ, params;
               restriction_ipeps = _restriction_ipeps)