function precondition_invese_single_envir(A, grad, rt, params, restriction_ipeps, fδEi, iter_precond)
    # size(A) == (1,) || throw(Base.error("precondition only supports 1x1 unit cell currently"))
    if fδEi[2] > 0.01 || fδEi[3] <= iter_precond
        return grad
    end
    δ = fδEi[2]
    A = restriction_ipeps(A)
    A = build_A(A, params)
    M = build_M(A, params) 
    # rt′ = leading_boundary(rt, M, params.boundary_alg)
    # Zygote.@ignore params.reuse_env && update!(rt, rt′)

    env = VUMPSEnv(rt, M, params.boundary_alg)
    @unpack ACu, ACd, FLo, FRo = env

    function re(x)
        s = size(x)
        χ, D = s[1], Int(sqrt(s[2]))
        reshape(x, χ, D, D, χ)
    end

    # @show size(grad[1])
    # @show size(ein"(((iaej,jbfk),abcdp),lcgk),idhl->efghp"(re(FLo[1]),re(conj(ACd[1])),grad[1],re(FRo[1]),re(ACu[1])))
    gradnew = deepcopy(grad)
    Ni = size(M)[1]
    forloop_iter = params.forloop_iter
    for p in 1:length(M)
        i, j = Tuple(findfirst(==(p), M.pattern))
        ir = Ni + 1 - i
        n = contract_n1(FLo[i,j], ACu[i,j], A[i,j], conj(ACd[ir,j]), FRo[i,j]; forloop_iter)
        go = (i + j) % 2 == 0 ? grad[:,:,:,:,:,p] : permutedims(grad[:,:,:,:,:,p], (3,4,1,2,5))
        # @show size(re(FLo[i,j])) size(re(conj(ACd[ir,j]))) size(re(FRo[i,j])) size(re(ACu[i,j]))
        if params.ifflatten 
            gnew, _ = linsolve(x->δ * x + ein"(((iaej,jbfk),abcdp),lcgk),idhl->efghp"(re(FLo[i,j]),re(conj(ACd[ir,j])),x,re(FRo[i,j]),re(ACu[i,j]))/n, go; isposdef = true, maxiter=1)
        else
            gnew, _ = linsolve(x->δ * x + ein"(((iaej,jbfk),abcdp),lcgk),idhl->efghp"(FLo[i,j],conj(ACd[ir,j]),x,FRo[i,j],ACu[i,j])/n, go; isposdef = true, maxiter=1)
        end
        gradnew[:,:,:,:,:,p] = (i + j) % 2 == 0 ? gnew : permutedims(gnew, (3,4,1,2,5))
    end

    return gradnew
end