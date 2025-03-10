function energy_value(model, A, M, env, params::iPEPSOptimize{:fermion, :square})
    atype = _arraytype(A[1])
    Ni, Nj = size(A)
    len = length(unique(A.pattern))
    D, d = space(A[1], 4), space(A[1], 3)
    DD = Zygote.@ignore fuse(D', D)
    SdD = Zygote.@ignore atype(swapgate(d, D))
    SDD = Zygote.@ignore atype(swapgate(D, D))
    h1, h2 = Zygote.@ignore atype.(hamiltonian(model))

    @unpack ACu, ARu, ACd, ARd, FLu, FRu, FLo, FRo = env

    It = Zygote.@ignore isomorphism(D'*D, DD)
    f1(x) = @tensor y[-1 -2 -3; -4] := It[-2 -3; 1] * x[-1 1; -4]
    f2(x) = @tensor y[-1 -2 -3; -4] := It'[1; -2 -3] * x[-1 1; -4]
    etol = 0
    loss = 0
    λ = params.λ
    ϵ = 0
    for p in 1:len
        i, j = Tuple(findfirst(==(p), A.pattern))
        params.verbosity >= 4 && println("===========$i,$j===========")
        ir = Ni + 1 - i
        jr = mod1(j + 1, Nj)
        e = horizontal_contraction(h1, h2, A[i,j],fdag(A[i,j], SDD),A[i,jr],fdag(A[i,jr], SDD),SdD,SDD,
                                    f1(FLo[i,j]),f1(ACu[i,j]),f1(ARu[i,jr]),f2(FRo[i,jr]),f1(ARd[ir,jr])',f1(ACd[ir,j])'
                                    )
        n = @tensor FLmap(FLmap(FLo[i,j], ACu[i,j], ACd[ir,j]', M[i,j]), ARu[i,jr], ARd[ir,jr]', M[i,jr])[1 2; 3] * FRo[i,jr][3 2; 1]
        params.verbosity >= 4 && println("h_H = $(e/n)")
        norm(n) < 1e-5 && @warn "n = $n"
        etol += e/n
        loss += e/(n + ϵ) + λ * 1/norm(n)

        ir  = mod1(i + 1, Ni)
        irr = mod1(Ni - i, Ni) 
        e = vertical_contraction(h1, h2, A[i,j],fdag(A[i,j], SDD),A[ir,j],fdag(A[ir,j], SDD),SdD,SDD,
                                  f1(ACu[i,j]),f2(FRu[i,j]),f2(FRo[ir,j]),f1(ACd[irr,j])',f1(FLo[ir,j]),f1(FLu[i,j])
                                 )
        n = @tensor ACmap(ACmap(ACu[i,j], FLu[i,j], FRu[i,j], M[i,j]), FLo[ir,j], FRo[ir,j], M[ir,j])[1 2; 3] * ACd[irr,j]'[3; 1 2]
        params.verbosity >= 4 && println("h_V = $(e/n)")
        norm(n) < 1e-5 && @warn "n = $n"
        etol += e/n
        loss += e/(n + ϵ) + λ * 1/norm(n)

        if params.ifNN
            ir  = mod1(i + 1,  Ni)
            irr = mod1(Ni - i, Ni) 
            jr  = mod1(j + 1,  Nj)
            e1 = next_neighbor1_contraction(h1, h2, A[i,j], fdag(A[i,j], SDD), A[i,jr], fdag(A[i,jr], SDD), A[ir,j], fdag(A[ir,j], SDD), A[ir,jr], fdag(A[ir,jr], SDD), SdD, SDD, f1(ACu[i,j]), f1(ARu[i,jr]), f1(FLu[i,j]), f1(FLo[ir,j]), f2(FRu[i,jr]), f2(FRo[ir,jr]), f1(ACd[irr,j])', f1(ARd[irr,jr])')
            e2 = next_neighbor2_contraction(h1, h2, A[i,j], fdag(A[i,j], SDD), A[i,jr], fdag(A[i,jr], SDD), A[ir,j], fdag(A[ir,j], SDD), A[ir,jr], fdag(A[ir,jr], SDD), SdD, SDD, f1(ACu[i,j]), f1(ARu[i,jr]), f1(FLu[i,j]), f1(FLo[ir,j]), f2(FRu[i,jr]), f2(FRo[ir,jr]), f1(ACd[irr,j])', f1(ARd[irr,jr])')
            n = four_site_contraction(M[i,j], M[i,jr], M[ir,j], M[ir,jr], ACu[i,j], ARu[i,jr], ACd[irr,j]', ARd[irr,jr]', FLu[i,j], FRu[i,jr], FLo[ir,j], FRo[ir,jr])
            norm(n) < 1e-5 && @warn "n = $n"
            params.verbosity >= 4 && println("h_NN1 = $(e1/n)\nh_NN2 = $(e2/n)")
            etol += (e1+e2)/n
        end
    end

    params.verbosity >= 3 && println("energy = $(etol/len)")
    return loss/len
end

function energy_value(model::Topological_Insulator, A, M, env, params::iPEPSOptimize{:fermion, :square})
    atype = _arraytype(A[1])
    Ni, Nj = size(A)
    Ni == Nj == 2 || throw(ArgumentError("Ni == Nj == 2 is required for Topological_Insulator"))
    len = length(unique(A.pattern))
    D, d = space(A[1], 4), space(A[1], 3)
    DD = Zygote.@ignore fuse(D', D)
    SdD = Zygote.@ignore atype(swapgate(d, D))
    SDD = Zygote.@ignore atype(swapgate(D, D))
    
    @unpack ACu, ARu, ACd, ARd, FLu, FRu, FLo, FRo = env
    
    It = Zygote.@ignore isomorphism(D'*D, DD)
    f1(x) = @tensor y[-1 -2 -3; -4] := It[-2 -3; 1] * x[-1 1; -4]
    f2(x) = @tensor y[-1 -2 -3; -4] := It'[1; -2 -3] * x[-1 1; -4]

    hN1  = Zygote.@ignore atype.(hamiltonian(Topological_Insulator(1.0, pi/4)))
    hN2  = Zygote.@ignore atype.(hamiltonian(Topological_Insulator(1.0, -pi/4)))
    hNN1 = Zygote.@ignore atype.(hamiltonian(Topological_Insulator(1/sqrt(2), 0)))
    hNN2 = Zygote.@ignore atype.(hamiltonian(Topological_Insulator(-1/sqrt(2), 0)))
    etol = 0
    loss = 0
    len = length(unique(A.pattern))
    λ = params.λ
    for p in 1:len
        i, j = Tuple(findfirst(==(p), A.pattern))
        if (i+j) % 2 == 0
            hh = hN1
            hv = hN2
            hNN = hNN1
        else
            hh = hN2
            hv = hN1
            hNN = hNN2
        end
        params.verbosity >= 4 && println("===========$i,$j===========")
        ir = Ni + 1 - i
        jr = mod1(j + 1, Nj)
        e = horizontal_contraction(hh[1], hh[2], A[i,j],fdag(A[i,j], SDD),A[i,jr],fdag(A[i,jr], SDD),SdD,SDD,
                                    f1(FLo[i,j]),f1(ACu[i,j]),f1(ARu[i,jr]),f2(FRo[i,jr]),f1(ARd[ir,jr])',f1(ACd[ir,j])'
                                    )
        n = @tensor FLmap(FLmap(FLo[i,j], ACu[i,j], ACd[ir,j]', M[i,j]), ARu[i,jr], ARd[ir,jr]', M[i,jr])[1 2; 3] * FRo[i,jr][3 2; 1]
        params.verbosity >= 4 && println("h_H = $(e/n)")
        norm(n) < 1e-5 && @warn "n = $n"
        etol += e/n
        loss += e/n + λ * 1/norm(n)

        ir  = mod1(i + 1, Ni)
        irr = mod1(Ni - i, Ni) 
        e = vertical_contraction(hv[1], hv[2], A[i,j],fdag(A[i,j], SDD),A[ir,j],fdag(A[ir,j], SDD),SdD,SDD,
                                  f1(ACu[i,j]),f2(FRu[i,j]),f2(FRo[ir,j]),f1(ACd[irr,j])',f1(FLo[ir,j]),f1(FLu[i,j])
                                 )
        n = @tensor ACmap(ACmap(ACu[i,j], FLu[i,j], FRu[i,j], M[i,j]), FLo[ir,j], FRo[ir,j], M[ir,j])[1 2; 3] * ACd[irr,j]'[3; 1 2]
        params.verbosity >= 4 && println("h_V = $(e/n)")
        norm(n) < 1e-5 && @warn "n = $n"
        etol += e/n
        loss += e/n + λ * 1/norm(n)

        ir  = mod1(i + 1, Ni)
        irr = mod1(Ni - i, Ni) 
        jr = mod1(j + 1, Nj)
        e1 = next_neighbor1_contraction(hNN[1], hNN[2], A[i,j], fdag(A[i,j], SDD), A[i,jr], fdag(A[i,jr], SDD), A[ir,j], fdag(A[ir,j], SDD), A[ir,jr], fdag(A[ir,jr], SDD), SdD, SDD, f1(ACu[i,j]), f1(ARu[i,jr]), f1(FLu[i,j]), f1(FLo[ir,j]), f2(FRu[i,jr]), f2(FRo[ir,jr]), f1(ACd[irr,j])', f1(ARd[irr,jr])')
        e2 = next_neighbor2_contraction(hNN[1], hNN[2], A[i,j], fdag(A[i,j], SDD), A[i,jr], fdag(A[i,jr], SDD), A[ir,j], fdag(A[ir,j], SDD), A[ir,jr], fdag(A[ir,jr], SDD), SdD, SDD, f1(ACu[i,j]), f1(ARu[i,jr]), f1(FLu[i,j]), f1(FLo[ir,j]), f2(FRu[i,jr]), f2(FRo[ir,jr]), f1(ACd[irr,j])', f1(ARd[irr,jr])')
        n = four_site_contraction(M[i,j], M[i,jr], M[ir,j], M[ir,jr], ACu[i,j], ARu[i,jr], ACd[irr,j]', ARd[irr,jr]', FLu[i,j], FRu[i,jr], FLo[ir,j], FRo[ir,jr])
        @show n
        norm(n) < 1e-5 && @warn "n = $n"
        params.verbosity >= 4 && println("h_NN1 = $(e1/n)\nh_NN2 = $(e2/n)")
        etol += (e1+e2)/n
        loss += (e1+e2)/n + λ * 1/norm(n)
    end

    params.verbosity >= 3 && println("energy = $(etol/len)")
    return loss/len
end