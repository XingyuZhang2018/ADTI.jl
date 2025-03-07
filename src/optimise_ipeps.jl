@kwdef mutable struct iPEPSOptimize{system, lattice}
    boundary_alg::VUMPS # boundary contraction algorithm
    reuse_env::Bool = Defaults.reuse_env # reuse the environment of the previous iteration
    verbosity::Int = Defaults.verbosity # verbosity level
    optimizer = Defaults.optimizer # optimization algorithm
    folder::String = Defaults.folder # folder to save the results
    show_every::Int = Defaults.show_every # show the results every `show_every` iterations
    save_every::Int = Defaults.save_every # save the results every `save_every` iterations
    ifflatten::Bool = true # flatten the center tensors before the optimization
    ifcheckpoint::Bool = false # use checkpointing
    ifNN::Bool = false # include the nearest-neighbour terms in the Hamiltonian
end

function optimise_ipeps(A::StructArray, model, χ::VectorSpace, params::iPEPSOptimize;
                        restriction_ipeps = _restriction_ipeps)
    D = space(A[1], 4)

    A′ = restriction_ipeps(A)
    M = build_M(A′, params)
    rt = VUMPSRuntime(M, χ, params.boundary_alg)
    function f(A) 
        A = restriction_ipeps(A)
        return params.ifcheckpoint ? real(checkpoint(energy, A, model, rt, params)) : real(energy(A, model, rt, params))
    end
    function fg(x)
        return f(x), StructArray(gradient(f, x)[1].data, A.pattern)
    end
    alg = params.optimizer
    t0 = time()
    x, f, g, numfg, normgradhistory = optimize(fg, A, alg; 
                                               inner = _inner, 
                                               finalize! = (x, f, g, iter)->_finalize!(x, f, g, iter, D, χ, params, t0)
    )
    return x, f, g, numfg, normgradhistory
end

_inner(x, dx1, dx2) = real(dot(dx1, dx2))
function _finalize!(x, f, g, iter, D, χ, params, t0)
    @unpack folder = params
    message = @sprintf("i = %5d\tt = %0.2f sec\tenergy = %.15f \tgnorm = %.3e\n", iter, time() - t0, f, norm(g))

    folder = rejoinpath(folder, "D$(D)_χ$(χ)")
    !(ispath(folder)) && mkpath(folder)
    if params.verbosity >= 3 && iter % params.show_every == 0
        printstyled(message; bold=true, color=:red)
        flush(stdout)

        logfile = open(rejoinpath(folder, "history.log"), "a")
        write(logfile, message)
        close(logfile)
    end
    if params.save_every != 0 && iter % params.save_every == 0
        save(rejoinpath(folder, "ipeps", "ipeps_No.$(iter).jld2"), "bcipeps", Array(x))
    end
    
    return x, f, g
end 