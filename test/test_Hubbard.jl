using ADIT
using CUDA
using LinearAlgebra
using Random
using TeneT
using Zygote
using U1ArrayKit

#####################################    parameters      ###################################
Random.seed!(100)
atype = Array
sitetype = electronZ2()
Ni, Nj = 1, 1
qnD = [0, 1]
dimsD = [1, 1]
qnχ = [0, 1]
dimsχ = [5, 5]

χ = sum(dimsχ)
No = 0
U1info = [qnD, qnχ, dimsD, dimsχ]
d = 4
model = Hubbard(1.0, 12.0, 6.0)
system = :fermion
lattice = :square
folder = "data/$model/$system/$lattice/$(Ni)x$(Nj)/"
boundary_alg = VUMPS(U1info=U1info,
                     ifupdown=false,
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
                                        sitetype=sitetype,
                                        reuse_env=true, 
                                        ifcheckpoint=true, 
                                        ifflatten=true,
                                        ifNN=false,
                                        verbosity=4, 
                                        maxiter=0,
                                        tol=1e-10,
                                        folder=folder
)
A = init_ipeps(;atype, No, d, Ni, Nj, params)
# A = ADIT.init_ipeps_spin111(;atype, model, params, No, ifWp=true, ϵ = 0, χ, Ni, Nj)
# A = ADIT.init_ipeps_h5(;atype, params, ifWp=false, ifreal=true, model, file="./data/kitsShf_sikh3nfr1D8D8.h5", D, Ni, Nj)
############################################################################################

function _restriction_ipeps(A)
    # Ar = Zygote.Buffer(A)
    # Ni, Nj = size(A)[[6,7]]
    # for j in 1:Nj, i in 1:Ni
    #     # Ar[:,:,:,:,:,i,j] = A[:,:,:,:,:,i,j] + permutedims(A[:,:,:,:,:,i,j], (2,1,3,5,4))
    #     if (i + j) % 2 == 0
    #         Ar[:,:,:,:,:,i,j] = A[:,:,:,:,:,1,1]
    #     else
    #         Ar[:,:,:,:,:,i,j] = A[:,:,:,:,:,1,2]
    #     end
    #     # Ar[:,:,:,:,:,i,j] = A[:,:,:,:,:,i,j] + permutedims(A[:,:,:,:,:,i,j], (2,1,3,5,4))
    # end
    # Ar = copy(Ar)
    # return Ar/norm(Ar)
    return A/norm(A)
end
optimise_ipeps(A, model, χ, params;
               restriction_ipeps = _restriction_ipeps)