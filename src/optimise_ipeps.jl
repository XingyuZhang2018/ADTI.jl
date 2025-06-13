@kwdef mutable struct iPEPSOptimize{system, lattice}
    pattern::Matrix{Int} # pattern of the iPEPS and the environment
    boundary_alg::VUMPS # boundary contraction algorithm
    reuse_env::Bool = Defaults.reuse_env # reuse the environment of the previous iteration
    
    optimizer = Defaults.optimizer # optimization algorithm
    ifprecondition::Bool = Defaults.ifprecondition
    iter_precond::Int = 20 # number of iterations for the preconditioner

    forloop_iter::Int = 1  # forloop_iter for energy calculation
    SUτ::Real = 0.0 # simple update τ 

    folder::String = Defaults.folder # folder to save the results
    show_every::Int = Defaults.show_every # show the results every `show_every` iterations
    save_every::Int = Defaults.save_every # save the results every `save_every` iterations
    ifsave_env::Bool = Defaults.ifsave_env # save the environment after the optimization
    ifload_env::Bool = Defaults.ifload_env # load the environment from the file

    ifflatten::Bool = true # flatten the center tensors before the optimization
    ifcheckpoint::Bool = false # use checkpointing
    ifNN::Bool = false # include the nearest-neighbour terms in the Hamiltonian

    verbosity::Int = Defaults.verbosity # verbosity level
end


"""
    optimise_ipeps(A::AbstractArray, model, χ::Int, params::iPEPSOptimize)

Optimise the iPEPS tensor `A` for the `model` Hamiltonian with bond dimension `χ`
using the parameters `params`. 
"""
function optimise_ipeps(A::AbstractArray, model, χ::Int, params::iPEPSOptimize;
                        restriction_ipeps = _restriction_ipeps)
    D = size(A, 1)
    rt = initialize_vumps_runtime(A, D, χ, params; restriction_ipeps)
    rt′ = deepcopy(rt)

    function f(A) 
        A = restriction_ipeps(A)
        return params.ifcheckpoint ? real(checkpoint(energy, A, model, rt, rt′, params)) : real(energy(A, model, rt, rt′, params))
    end
    function fg(x)
        # return f(x), gradient(f, x)[1]
        e, vjp = pullback(f, x)
        return e, vjp(1)[1]
    end
    alg = params.optimizer
    t0 = time()
    fδEi = [1.0,1.0,0]
    _precondition(x, g) = params.ifprecondition ? precondition_invese_single_envir(x, g, rt, params, restriction_ipeps, fδEi, params.iter_precond) : g
    x, f, g, numfg, normgradhistory = optimize(fg, A, alg; 
                                               precondition=_precondition, 
                                               inner = _inner, 
                                               finalize! = (x, f, g, iter)->_finalize!(x, f, g, iter, rt, rt′, D, χ, params, t0, fδEi)
    )
    return x, f, g, numfg, normgradhistory
end

_inner(x, dx1, dx2) = real(dot(dx1, dx2))
function _finalize!(x, f, g, iter, rt, rt′, D, χ, params, t0, fδEi)
    @unpack folder = params

    fδEi[3] = iter
    fδEi[2] = abs(fδEi[1] - f)
    fδEi[1] = f
    message = @sprintf("i = %5d\tt = %0.2f sec\tenergy = %.15f \tgnorm = %.3e\n", iter, time() - t0, f, norm(g))

    folder = joinpath(folder, "D$(D)_χ$(χ)")
    !(ispath(folder)) && mkpath(folder)
    params.reuse_env && update!(rt, rt′)
    params.ifsave_env && save_rt(folder, rt)

    if params.verbosity >= 3 && iter % params.show_every == 0
        printstyled(message; bold=true, color=:red)
        flush(stdout)

        logfile = open(joinpath(folder, "history.log"), "a")
        write(logfile, message)
        close(logfile)
    end
    if params.save_every != 0 && iter % params.save_every == 0
        save(joinpath(folder, "ipeps", "ipeps_No.$(iter).jld2"), "bcipeps", Array(x))
    end
    
    return x, f, g
end 