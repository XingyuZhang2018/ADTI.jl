function energy_value(model, A, env, params::iPEPSOptimize{:boson, :merge})
    @unpack ACu, ARu, ACd, ARd, FLu, FRu, FLo, FRo = env
    atype = _arraytype(ACu[1])
    S = model.S
    d = Int(2*S + 1)
    Sx = Zygote.@ignore const_Sx(S)
    Sy = Zygote.@ignore const_Sy(S)
    Sz = Zygote.@ignore const_Sz(S)

    Ni, Nj = size(ACu)
    etol = 0
    for j = 1:Nj, i = 1:Ni
        if (i,j) in [(1,1),(2,2)]
            Jx, Jy, Jz = model.Jx*1.0, model.Jy*1.0, model.Jz
        elseif (i,j) in [(1,2),(2,3)]
            Jx, Jy, Jz = model.Jx*1.0, model.Jy, model.Jz*1.0
        else
            Jx, Jy, Jz = model.Jx, model.Jy*1.0, model.Jz*1.0
        end
        params.verbosity >= 4 && println("===========$i,$j===========")
        ir = Ni + 1 - i
        jr = mod1(j + 1, Nj)
        O1 = atype(Jx * reshape(ein"ac,bd->abcd"(I(d), Sx), d^2,d^2))
        O2 = atype(reshape(ein"ac,bd->abcd"(Sx, I(d)), d^2,d^2))
        e = contract_o2_H(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,jr],ARu[i,jr],A[i,jr],conj(ARd[ir,jr]), O1, O2)
        n = contract_n2_H(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,jr],ARu[i,jr],A[i,jr],conj(ARd[ir,jr]))
        params.verbosity >= 4 && println("hx = $(e/n)")
        etol += e/n

        O = atype(Jx * reshape(ein"ac,bd->abcd"(Sz, Sz), d^2,d^2))
        e = contract_o1(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,j], O)
        n = contract_n1(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,j])
        etol += e/n
        params.verbosity >= 4 && println("hz = $(e/n)")
        
        ir  = mod1(i + 1, Ni)
        irr = mod1(Ni - i, Ni) 
        O1 = atype(Jy * reshape(ein"ac,bd->abcd"(I(d), Sy), d^2,d^2))
        O2 = atype(reshape(ein"ac,bd->abcd"(Sy, I(d)), d^2,d^2))
        e = contract_o2_V(ACu[i,j],FLu[i,j],A[i,j],FRu[i,j],FLo[ir,j],A[ir,j],FRo[ir,j],conj(ACd[irr,j]), O1, O2)
        n = contract_n2_V(ACu[i,j],FLu[i,j],A[i,j],FRu[i,j],FLo[ir,j],A[ir,j],FRo[ir,j],conj(ACd[irr,j]))
        params.verbosity >= 4 && println("hy = $(e/n)")
        etol += e/n
    end

    params.verbosity >= 3 && println("energy = $(etol/Ni/Nj/2)")
    return etol/Ni/Nj/2
end

function energy_value(model, A, env, params::iPEPSOptimize{:boson, :brickwall})
    @unpack ACu, ARu, ACd, ARd, FLu, FRu, FLo, FRo = env
    atype = _arraytype(ACu[1])
    S = model.S
    Sx = Zygote.@ignore const_Sx(S)
    Sy = Zygote.@ignore const_Sy(S)
    Sz = Zygote.@ignore const_Sz(S)
    Sz2 = Sz * Sz 

    Ni, Nj = size(ACu)
    etol = 0
    for j = 1:Nj, i = 1:Ni
        params.verbosity >= 4 && println("===========$i,$j===========")
        if (i + j) % 2 != 0
            ir  = mod1(i + 1, Ni)
            irr = mod1(Ni - i, Ni) 
            e = contract_o2_V(ACu[i,j],FLu[i,j],A[i,j],FRu[i,j],FLo[ir,j],A[ir,j],FRo[ir,j],conj(ACd[irr,j]), atype(model.Jy * Sy), atype(Sy))
            n = contract_n2_V(ACu[i,j],FLu[i,j],A[i,j],FRu[i,j],FLo[ir,j],A[ir,j],FRo[ir,j],conj(ACd[irr,j]))
            params.verbosity >= 4 && println("hy = $(e/n)")
            etol += e/n

            O_H = sqrt(model.Jx) * Sx
        else
            O_H = sqrt(model.Jz) * Sz
        end

        ir = Ni + 1 - i
        jr = mod1(j + 1, Nj)
        e = contract_o2_H(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,jr],ARu[i,jr],A[i,jr],conj(ARd[ir,jr]), atype(O_H), atype(O_H))
        n = contract_n2_H(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,jr],ARu[i,jr],A[i,jr],conj(ARd[ir,jr]))
        params.verbosity >= 4 && (i + j) % 2 != 0 ? println("hz = $(e/n)") : println("hx = $(e/n)")
        etol += e/n
    end

    params.verbosity >= 3 && println("energy = $(etol/Ni/Nj)")
    return etol/Ni/Nj
end

function magnetization_value(model, A, env, params::iPEPSOptimize{:boson, :merge})
    @unpack ACu, ARu, ACd, ARd, FLu, FRu, FLo, FRo = env
    atype = _arraytype(ACu[1])
    S = model.S
    d = Int(2*S + 1)
    Sx = Zygote.@ignore const_Sx(S)
    Sy = Zygote.@ignore const_Sy(S)
    Sz = Zygote.@ignore const_Sz(S)

    Ni, Nj = size(ACu)
    Mag = Array{Array{ComplexF64,1},3}(undef, Ni, Nj, 2)
    Mnorm = Array{ComplexF64,3}(undef, Ni, Nj, 2)
    for j = 1:Nj, i = 1:Ni
        params.verbosity >= 4 && println("===========$i,$j===========")
        ir = Ni + 1 - i
        O1 = atype(reshape(ein"ac,bd->abcd"(Sx, I(d)), d^2,d^2))
        O2 = atype(reshape(ein"ac,bd->abcd"(I(d), Sx), d^2,d^2))
        Mx1 = contract_o1(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,j], O1)
        Mx2 = contract_o1(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,j], O2)

        O1 = atype(reshape(ein"ac,bd->abcd"(Sy, I(d)), d^2,d^2))
        O2 = atype(reshape(ein"ac,bd->abcd"(I(d), Sy), d^2,d^2))
        My1 = contract_o1(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,j], O1)
        My2 = contract_o1(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,j], O2)

        O1 = atype(reshape(ein"ac,bd->abcd"(Sz, I(d)), d^2,d^2))
        O2 = atype(reshape(ein"ac,bd->abcd"(I(d), Sz), d^2,d^2))
        Mz1 = contract_o1(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,j], O1)
        Mz2 = contract_o1(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,j], O2)

        n = contract_n1(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,j])
        Mag[i,j,1] = [Mx1/n, My1/n, Mz1/n]
        Mag[i,j,2] = [Mx2/n, My2/n, Mz2/n]
        Mnorm[i,j,1] = norm(Mag[i,j,1])
        Mnorm[i,j,2] = norm(Mag[i,j,2])
        params.verbosity >= 4 && println("M1 = $(Mag[i,j,1])\nM2 = $(Mag[i,j,2])\n|M1| = $(Mnorm[i,j,1])\n|M2| = $(Mnorm[i,j,2])")
    end
    params.verbosity >= 4 && println("|M|_mean = $(sum(Mnorm)/Ni/Nj/2)")

    return Mag, Mnorm
end

function magnetization_value(model, A, env, params::iPEPSOptimize{:boson, :brickwall})
    @unpack ACu, ARu, ACd, ARd, FLu, FRu, FLo, FRo = env
    atype = _arraytype(ACu[1])
    S = model.S
    Sx = Zygote.@ignore const_Sx(S)
    Sy = Zygote.@ignore const_Sy(S)
    Sz = Zygote.@ignore const_Sz(S)

    Ni, Nj = size(ACu)
    Mag = Array{Array{ComplexF64,1},2}(undef, Ni, Nj)
    Mnorm = Array{ComplexF64,2}(undef, Ni, Nj)
    for j = 1:Nj, i = 1:Ni
        params.verbosity >= 4 && println("===========$i,$j===========")
        ir = Ni + 1 - i
        Mx = contract_o1(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,j], atype(Sx))
        My = contract_o1(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,j], atype(Sy))
        Mz = contract_o1(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,j], atype(Sz))
        
        n = contract_n1(FLo[i,j],ACu[i,j],A[i,j],conj(ACd[ir,j]),FRo[i,j])
        Mag[i,j] = [Mx/n, My/n, Mz/n]
        Mnorm[i,j] = norm(Mag[i,j])
        params.verbosity >= 4 && println("M = $(Mag[i,j])\n|M| = $(Mnorm[i,j])")
    end

    params.verbosity >= 4 && println("|M|_mean = $(sum(Mnorm)/Ni/Nj)")
    return Mag, Mnorm
end

