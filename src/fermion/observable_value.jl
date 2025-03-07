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
        
        ir  = mod1(i + 1, Ni)
        irr = mod1(Ni - i, Ni) 
        e = vertical_contraction(h1, h2, A[i,j],fdag(A[i,j], SDD),A[ir,j],fdag(A[ir,j], SDD),SdD,SDD,
                                  f1(ACu[i,j]),f2(FRu[i,j]),f2(FRo[ir,j]),f1(ACd[irr,j])',f1(FLo[ir,j]),f1(FLu[i,j])
                                 )
        n = @tensor ACmap(ACmap(ACu[i,j], FLu[i,j], FRu[i,j], M[i,j]), FLo[ir,j], FRo[ir,j], M[ir,j])[1 2; 3] * ACd[irr,j]'[3; 1 2]
        params.verbosity >= 4 && println("h_V = $(e/n)")
        norm(n) < 1e-5 && @warn "n = $n"
        etol += e/n
    end

    params.verbosity >= 3 && println("energy = $(etol/len)")
    return etol/len
end

function energy_value(model::Topological_Insulator, A, env, params::iPEPSOptimize{:fermion, :square})
    A, Adag = A
    atype = _arraytype(A[1])
    Ni, Nj = size(A)
    Ni == Nj == 2 || throw(ArgumentError("Ni == Nj == 2 is required for Topological_Insulator"))
    D, d = size(A[1])[[2,3]]
    SdD = atype(swapgate(d, D))
    SDD = atype(swapgate(D, D))
    @unpack ACu, ARu, ACd, ARd, FLu, FRu, FLo, FRo = env

    hN1  = atype.(hamiltonian(Topological_Insulator(1.0, pi/4)))
    hN2  = atype.(hamiltonian(Topological_Insulator(1.0, -pi/4)))
    hNN1 = atype.(hamiltonian(Topological_Insulator(1/sqrt(2), 0)))
    hNN2 = atype.(hamiltonian(Topological_Insulator(-1/sqrt(2), 0)))
    etol = 0
    for j = 1:Nj, i = 1:1
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
        ρ = oc_H_fermion(A[i,j],SdD,Adag[i,j],SdD,SDD,SDD,Adag[i,jr],SdD,SdD,SDD,A[i,jr],SDD,
                         FLo[i,j],ACu[i,j],ARu[i,jr],FRo[i,jr],conj(ARd[ir,jr]),conj(ACd[ir,j])
        )
        e = sum(ein"ijkl,ijkl -> "(ρ,hh))
        n = sum(ein"ijij->"(ρ))
        params.verbosity >= 4 && println("h_H = $(e/n)")
        etol += e/n
        @show n
        ir  = mod1(i + 1, Ni)
        irr = mod1(Ni - i, Ni) 
        ρ = oc_V_fermion(A[i,j],Adag[i,j],SDD,SdD,SdD,SDD,A[ir,j],Adag[ir,j],SDD,SdD,SdD,SDD,
                         ACu[i,j],FRu[i,j],FRo[ir,j],conj(ACd[irr,j]),FLo[ir,j],FLu[i,j]
        )
        e = sum(ein"ijkl,ijkl -> "(ρ,hv))
        n = sum(ein"ijij->"(ρ))
        params.verbosity >= 4 && println("h_V = $(e/n)")
        etol += e/n
        @show n
        ir  = mod1(i + 1, Ni)
        irr = mod1(Ni - i, Ni) 
        jr = mod1(j + 1, Nj)
        ρ = oc_NN1_fermion(SDD,A[i,j],SdD,SdD,SDD,SdD,SdD,Adag[i,j],
                        SDD,A[i,jr],SDD,Adag[i,jr],
                        SDD,A[ir,j],SDD,Adag[ir,j],
                        SdD,SdD,SDD,SdD,SdD,A[ir,jr],SDD,Adag[ir,jr],
                        ACu[i,j],ARu[i,jr],FLu[i,j],FLo[ir,j],FRu[i,jr],FRo[ir,jr],conj(ACd[irr,j]),conj(ARd[irr,jr])
        )
        e = sum(ein"ijkl,ijkl -> "(ρ,hNN))
        n = sum(ein"ijij->"(ρ))
        params.verbosity >= 4 && println("h_NN1 = $(e/n)")
        etol += e/n
        @show n
        ρ = oc_NN2_fermion(SDD,A[i,j],SDD,Adag[i,j],
                        SDD,A[i,jr],SDD,Adag[i,jr],SdD,SdD,SdD,SdD,
                        SDD,A[ir,j],SdD,SdD,SdD,SdD,SDD,Adag[ir,j],
                        SDD,A[ir,jr],SDD,Adag[ir,jr],
                        ACu[i,j],ARu[i,jr],FLu[i,j],FLo[ir,j],FRu[i,jr],FRo[ir,jr],conj(ACd[irr,j]),conj(ARd[irr,jr])
        )
        e = sum(ein"ijkl,ijkl -> "(ρ,hNN))
        n = sum(ein"ijij->"(ρ))
        params.verbosity >= 4 && println("h_NN2 = $(e/n)")
        etol += e/n
        @show n
    end

    params.verbosity >= 3 && println("energy = $(etol/Ni/Nj*2)")
    return etol/Ni/Nj*2
end