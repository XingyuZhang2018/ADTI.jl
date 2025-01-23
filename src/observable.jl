"""
    energy(h, bcipeps; χ, tol, maxiter)
return the energy of the `bcipeps` 2-site hamiltonian `h` and calculated via a
BCVUMPS with parameters `χ`, `tol` and `maxiter`.
"""
function energy(A, model, rt, params::iPEPSOptimize)
    A = bulid_A(A, params)
    M = bulid_M(A, params)
    rt′ = leading_boundary(rt, M, params.boundary_alg)
    Zygote.@ignore params.reuse_env && update!(rt, rt′)
    env = VUMPSEnv(rt′, M, params.boundary_alg)
    return energy_value(model, A, env, params)
end

function observable(A, model, χ, params::iPEPSOptimize;
                    restriction_ipeps = _restriction_ipeps)
    A = restriction_ipeps(A)
    A = bulid_A(A, params)
    M = bulid_M(A, params)
    rt = VUMPSRuntime(M, χ, params.boundary_alg)
    rt = leading_boundary(rt, M, params.boundary_alg)
    env = VUMPSEnv(rt, M, params.boundary_alg)

    e = energy_value(model, A, env, params)
    mag = magnetization_value(model, A, env, params)

    return e, mag
end
