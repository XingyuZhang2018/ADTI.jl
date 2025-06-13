function initialize_vumps_runtime(A, D, χ, params; restriction_ipeps)
    # Construct the expected file path
    folder_path = joinpath(params.folder, "D$(D)_χ$(χ)")
    file_path = joinpath(folder_path, "VUMPS_rt_env.jld2")

    if params.ifload_env
        if ispath(file_path)
            try
                return load_rt(folder_path, _arraytype(A))
            catch e
                @warn "Failed to load runtime environment from $file_path: $(sprint(showerror, e)). Creating new environment."
                return create_new_runtime(A, χ, params; restriction_ipeps)
            end
        else
            @warn "File $file_path does not exist. Creating new VUMPS runtime environment."
            return create_new_runtime(A, χ, params; restriction_ipeps)
        end
    else
        return create_new_runtime(A, χ, params; restriction_ipeps)
    end
end

function create_new_runtime(A, χ, params; restriction_ipeps)
    A = restriction_ipeps(A)
    A = build_A(A, params)
    M = build_M(A, params)
    return VUMPSRuntime(M, χ, params.boundary_alg)
end
