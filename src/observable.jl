"""
    energy(h, bcipeps; χ, tol, maxiter)
return the energy of the `bcipeps` 2-site hamiltonian `h` and calculated via a
BCVUMPS with parameters `χ`, `tol` and `maxiter`.
"""
function energy(A, model, rt, rt′, params::iPEPSOptimize)
    A = build_A(A, params)
    M = build_M(A, params)
    rt, _ = leading_boundary(rt, M, params.boundary_alg)
    Zygote.@ignore update!(rt′, rt)
    env = VUMPSEnv(rt, M, params.boundary_alg)
    return energy_value(model, A, M, env, params)
end

function observable(A, model, χ, params::iPEPSOptimize;
                    restriction_ipeps = _restriction_ipeps)
    rt = initialize_vumps_runtime(A, model, χ, params; restriction_ipeps)
    rt = leading_boundary(rt, M, params.boundary_alg)
    env = VUMPSEnv(rt, M, params.boundary_alg)

    e = energy_value(model, A, env, params)
    mag = magnetization_value(model, A, env, params)

    return e, mag
end
