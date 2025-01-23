@kwdef mutable struct iPEPSOptimize{system, lattice}
    boundary_alg::VUMPS # boundary contraction algorithm
    reuse_env::Bool = Defaults.reuse_env # reuse the environment of the previous iteration
    verbosity::Int = Defaults.verbosity # verbosity level
    maxiter::Int = Defaults.fpgrad_maxiter # maximum number of optimization iterations
    tol::Real = Defaults.fpgrad_tol # tolerance for the optimization
    optimizer = Defaults.optimizer # optimization algorithm
    folder::String = Defaults.folder # folder to save the results
    show_every::Int = Defaults.show_every # show the results every `show_every` iterations
    save_every::Int = Defaults.save_every # save the results every `save_every` iterations
    ifflatten::Bool = true # flatten the center tensors before the optimization
    ifcheckpoint::Bool = false # use checkpointing
    ifNN::Bool = false # include the nearest-neighbour terms in the Hamiltonian
end


"""
    optimise_ipeps(A::AbstractArray, model, χ::Int, params::iPEPSOptimize)

Optimise the iPEPS tensor `A` for the `model` Hamiltonian with bond dimension `χ`
using the parameters `params`. 
"""
function optimise_ipeps(A::AbstractArray, model, χ::Int, params::iPEPSOptimize;
                        restriction_ipeps = _restriction_ipeps)
    D = size(A, 1)
    A′ = restriction_ipeps(A)
    A′ = bulid_A(A′, params)
    M = bulid_M(A′, params)
    rt = VUMPSRuntime(M, χ, params.boundary_alg)

    function f(A) 
        A = restriction_ipeps(A)
        return params.ifcheckpoint ? real(checkpoint(energy, A, model, rt, params)) : real(energy(A, model, rt, params))
    end
    function g(A)
        # f(x)
        grad = Zygote.gradient(f,A)[1]
        return grad
    end
    res = optimize(f, g, 
        A, params.optimizer, inplace = false,
        Optim.Options(f_tol=params.tol, 
                      iterations=params.maxiter,
                      extended_trace=true,
                      callback=os->writelog(os, params, D, χ)
        )
    )
    return res
end

"""
    writelog(os::OptimizationState, key=nothing)

Write the current state of the optimization `os` to a log file. The log file is
saved in the folder `key` with the name `history.log`. The log file contains the
iteration number, the time spent on the iteration, the energy and the gradient norm.
The log file is appended to if it already exists.

"""
function writelog(os::OptimizationState, params::iPEPSOptimize, D::Int, χ::Int)
    @unpack folder = params

    message = @sprintf("i = %5d\tt = %0.2f sec\tenergy = %.15f \tgnorm = %.3e\n", os.iteration, os.metadata["time"], os.value, os.g_norm)

    folder = joinpath(folder, "D$(D)_χ$(χ)")
    !(ispath(folder)) && mkpath(folder)
    if params.verbosity >= 3 && os.iteration % params.show_every == 0
        printstyled(message; bold=true, color=:red)
        flush(stdout)

        logfile = open(joinpath(folder, "history.log"), "a")
        write(logfile, message)
        close(logfile)
    end
    if params.save_every != 0 && os.iteration % params.save_every == 0
        save(joinpath(folder, "ipeps", "ipeps_No.$(os.iteration).jld2"), "bcipeps", Array(os.metadata["x"]))
    end

    return false
end