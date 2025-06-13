@kwdef mutable struct iPEPSOptimize{system, lattice}
    boundary_alg::VUMPS # boundary contraction algorithm
    reuse_env::Bool = Defaults.reuse_env # reuse the environment of the previous iteration
    λ::Real = 1e-3 # regularization parameter for not 0/0 energy
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
    function fn(x)
        A = restriction_ipeps(x)
        return overlap(A, rt, params)
        # M = build_M(A, params)
        # loss = norm(M.data[1] - _fit_spaces(permute(adjoint(M.data[1]), ((1,4), (3,2))), M.data[1])) +  
        #         norm(M.data[1] - _fit_spaces(permute(adjoint(M.data[1]), ((3,2), (1,4))), M.data[1]))
        # @show loss
        # return loss
    end
    function fg(x)
        g = StructArray(gradient(f, x)[1].data, A.pattern)
        gn = StructArray(gradient(fn, x)[1].data, A.pattern)
        α = dot(g, gn) / dot(gn, gn)
        @show dot(g, gn) / norm(gn)/ norm(g)
        return f(x), g 
    end
    alg = params.optimizer
    t0 = time()
    e_diff = [0.0, Inf]
    # function precondition(x, g) 
    #     gn = StructArray(gradient(fn, x)[1].data, A.pattern)
    #     g_new = g - dot(g, gn) / dot(gn, gn) * gn
    #     return g_new
    # end
    x, f, g, numfg, normgradhistory = optimize(fg, A, alg; 
                                            #    precondition, 
                                               inner = _inner, 
                                               finalize! = (x, f, g, iter)->_finalize!(x, f, g, iter, D, χ, params, t0, e_diff)
    )
    return x, f, g, numfg, normgradhistory
end

_inner(x, dx1, dx2) = real(dot(dx1, dx2))
function _finalize!(x, f, g, iter, D, χ, params, t0, e_diff)
    @unpack folder = params
    # M = build_M(x/norm(x), params)
    # env = VUMPSEnv(rt, M, params.boundary_alg)
    # @unpack ACu, ACd, FLo, FRo = env
    # n = @tensor FLmap(FLo[1], ACu[1], ACd[1]', M[1])[1 2; 3] * FRo[1][3 2; 1]
    # message = @sprintf("i = %5d\tt = %0.2f sec\tenergy = %.15f \tgnorm = %.3e\tλ = %.3e\tn = %.15f\n", iter, time() - t0, f, norm(g), params.λ, norm(n))
    message = @sprintf("i = %5d\tt = %0.2f sec\tenergy = %.15f \tgnorm = %.3e\tλ = %.3e\n", iter, time() - t0, f, norm(g), params.λ)
    # e_diff[2] = norm(e_diff[1] - f)
    # e_diff[1] = f
    # if e_diff[2] < 1e-4 &&  params.λ > 1e-5
    #     params.λ /= 10
    # end
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