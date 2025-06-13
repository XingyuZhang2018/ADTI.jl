module ADTI
using BitBasis
using ChainRulesCore
using FileIO
using KrylovKit
using LinearAlgebra
using LineSearches
using Random
using TensorKit
using OptimKit
using OMEinsum
using Printf
using Parameters
using Zygote
using TeneT
using JLD2

using TeneT: _arraytype, update!, checkpoint, leg3, leg4, to_CuArray
using TeneT: ALCtoAC, FLmap, ACmap
using TeneT: _fit_spaces

export Kitaev, Hubbard, Free_Fermion, Topological_Insulator
export observable
export iPEPSOptimize, init_ipeps, optimise_ipeps, energy

include("defaults.jl")
include("patch.jl")
include("hamiltonian_models.jl")
include("optimise_ipeps.jl")
include("init_ipeps.jl")
# include("boson/build_A_M.jl")
# include("boson/contraction.jl")
# include("boson/observable_value.jl")
include("fermion/build_M.jl")
include("fermion/contraction.jl")
include("fermion/observable_value.jl")
include("observable.jl")

end
