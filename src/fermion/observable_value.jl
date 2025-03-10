function energy_value(model, A, env, params::iPEPSOptimize{:fermion, :square})
    A, Adag = A
    atype = _arraytype(A[1])
    Ni, Nj = size(A)
    D, d = size(A[1])[[2,3]]
    SdD = atype(swapgate(d, D))
    SDD = atype(swapgate(D, D))
    h = atype(hamiltonian(model))

    @unpack ACu, ARu, ACd, ARd, FLu, FRu, FLo, FRo = env
    # χ = size(ACu[1], 1)
    # ACu, ARu, ACd, ARd, FLu, FRu, FLo, FRo = map(x->map(y->reshape(y, χ,D,D,χ), x), [ACu, ARu, ACd, ARd, FLu, FRu, FLo, FRo])

    etol = 0
    for j = 1:Nj, i = 1:Ni
        params.verbosity >= 4 && println("===========$i,$j===========")
        ir = Ni + 1 - i
        jr = mod1(j + 1, Nj)
        ρ = oc_H_fermion(A[i,j],SdD,Adag[i,j],SdD,SDD,SDD,Adag[i,jr],SdD,SdD,SDD,A[i,jr],SDD,
                         FLo[i,j],ACu[i,j],ARu[i,jr],FRo[i,jr],conj(ARd[ir,jr]),conj(ACd[ir,j])
        )
        e = sum(ein"ijkl,ijkl -> "(ρ,h))
        n = sum(ein"ijij->"(ρ))
        params.verbosity >= 4 && println("h_H = $(e/n)")
        etol += e/n
        
        ir  = mod1(i + 1, Ni)
        irr = mod1(Ni - i, Ni) 
        ρ = oc_V_fermion(A[i,j],Adag[i,j],SDD,SdD,SdD,SDD,A[ir,j],Adag[ir,j],SDD,SdD,SdD,SDD,
                         ACu[i,j],FRu[i,j],FRo[ir,j],conj(ACd[irr,j]),FLo[ir,j],FLu[i,j]
        )
        e = sum(ein"ijkl,ijkl -> "(ρ,h))
        n = sum(ein"ijij->"(ρ))
        params.verbosity >= 4 && println("h_V = $(e/n)")
        etol += e/n

        if params.ifNN
            ir  = mod1(i + 1, Ni)
            irr = mod1(Ni - i, Ni) 
            jr = mod1(j + 1, Nj)
            ρ = oc_NN1_fermion(SDD,A[i,j],SdD,SdD,SDD,SdD,SdD,Adag[i,j],
                            SDD,A[i,jr],SDD,Adag[i,jr],
                            SDD,A[ir,j],SDD,Adag[ir,j],
                            SdD,SdD,SDD,SdD,SdD,A[ir,jr],SDD,Adag[ir,jr],
                            ACu[i,j],ARu[i,jr],FLu[i,j],FLo[ir,j],FRu[i,jr],FRo[ir,jr],conj(ACd[irr,j]),conj(ARd[irr,jr])
            )
            e = sum(ein"ijkl,ijkl -> "(ρ,h))
            n = sum(ein"ijij->"(ρ))
            params.verbosity >= 4 && println("h_NN1 = $(e/n)")
            etol += e/n

            ρ = oc_NN2_fermion(SDD,A[i,j],SDD,Adag[i,j],
                            SDD,A[i,jr],SDD,Adag[i,jr],SdD,SdD,SdD,SdD,
                            SDD,A[ir,j],SdD,SdD,SdD,SdD,SDD,Adag[ir,j],
                            SDD,A[ir,jr],SDD,Adag[ir,jr],
                            ACu[i,j],ARu[i,jr],FLu[i,j],FLo[ir,j],FRu[i,jr],FRo[ir,jr],conj(ACd[irr,j]),conj(ARd[irr,jr])
            )
            e = sum(ein"ijkl,ijkl -> "(ρ,h))
            n = sum(ein"ijij->"(ρ))
            params.verbosity >= 4 && println("h_NN2 = $(e/n)")
            etol += e/n
        end

    end

    params.verbosity >= 3 && println("energy = $(etol/Ni/Nj)")
    return etol/Ni/Nj
end

function energy_value(model::Topological_Insulator, A, M, env, params::iPEPSOptimize{:fermion, :square})
    A, Adag = A
    atype = _arraytype(A[1])
    Ni, Nj = size(A)
    Ni == Nj == 2 || throw(ArgumentError("Ni == Nj == 2 is required for Topological_Insulator"))
    D, d = size(A[1])[[2,3]]
    SdD = atype(swapgate(d, D))
    SDD = atype(swapgate(D, D))
    @unpack ACu, ARu, ACd, ARd, FLu, FRu, FLo, FRo = env

    hN1  = Zygote.@ignore atype.(hamiltonian(Topological_Insulator(1.0, pi/4)))
    hN2  = Zygote.@ignore atype.(hamiltonian(Topological_Insulator(1.0, -pi/4)))
    hNN1 = Zygote.@ignore atype.(hamiltonian(Topological_Insulator(1/sqrt(2), 0)))
    hNN2 = Zygote.@ignore atype.(hamiltonian(Topological_Insulator(-1/sqrt(2), 0)))
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
        e = sum(oc_H_fermion(A[i,j],SdD,Adag[i,j],SdD,SDD,SDD,Adag[i,jr],SdD,SdD,SDD,A[i,jr],SDD,hh[1],hh[2],
                         FLo[i,j],ACu[i,j],ARu[i,jr],FRo[i,jr],conj(ARd[ir,jr]),conj(ACd[ir,j])
        ))
        n = sum(ein"abcd,abcd->"(FLmap(FLmap(FLo[i,j], ACu[i,j], conj(ACd[ir,j]), M[i,j]), ARu[i,jr], conj(ARd[ir,jr]), M[i,jr]), FRo[i,jr]))
        params.verbosity >= 4 && println("h_H = $(e/n)")
        etol += e/n
        
        ir  = mod1(i + 1, Ni)
        irr = mod1(Ni - i, Ni) 
        e = sum(oc_V_fermion(A[i,j],Adag[i,j],SDD,SdD,SdD,SDD,A[ir,j],Adag[ir,j],SDD,SdD,SdD,SDD,hv[1],hv[2],
                         ACu[i,j],FRu[i,j],FRo[ir,j],conj(ACd[irr,j]),FLo[ir,j],FLu[i,j]
        ))
        n = sum(ein"abcd,abcd->"(ACmap(ACmap(ACu[i,j], FLu[i,j], FRu[i,j], M[i,j]), FLo[ir,j], FRo[ir,j], M[ir,j]), conj(ACd[irr,j])))
        params.verbosity >= 4 && println("h_V = $(e/n)")
        etol += e/n

        ir  = mod1(i + 1, Ni)
        irr = mod1(Ni - i, Ni) 
        jr = mod1(j + 1, Nj)
        e = sum(oc_NN1_fermion(SDD,A[i,j],SdD,SdD,SDD,SdD,SdD,Adag[i,j],
                        SDD,A[i,jr],SDD,Adag[i,jr],
                        SDD,A[ir,j],SDD,Adag[ir,j],
                        SdD,SdD,SDD,SdD,SdD,A[ir,jr],SDD,Adag[ir,jr],hNN[1],hNN[2],
                        ACu[i,j],ARu[i,jr],FLu[i,j],FLo[ir,j],FRu[i,jr],FRo[ir,jr],conj(ACd[irr,j]),conj(ARd[irr,jr])
        ))
        n = sum(oc_NN_n_fermion(FLu[i,j],ACu[i,j],M[i,j],ARu[i,jr],FRu[i,jr],M[i,jr],FLo[ir,j],conj(ACd[irr,j]),M[ir,j],conj(ARd[irr,jr]),FRo[ir,jr],M[ir,jr]))
        params.verbosity >= 4 && println("h_NN1 = $(e/n)")
        etol += e/n

        e = sum(oc_NN2_fermion(SDD,A[i,j],SDD,Adag[i,j],
                        SDD,A[i,jr],SDD,Adag[i,jr],SdD,SdD,SdD,SdD,
                        SDD,A[ir,j],SdD,SdD,SdD,SdD,SDD,Adag[ir,j],
                        SDD,A[ir,jr],SDD,Adag[ir,jr],hNN[1],hNN[2],
                        ACu[i,j],ARu[i,jr],FLu[i,j],FLo[ir,j],FRu[i,jr],FRo[ir,jr],conj(ACd[irr,j]),conj(ARd[irr,jr])
        ))
        params.verbosity >= 4 && println("h_NN2 = $(e/n)")
        etol += e/n
    end

    params.verbosity >= 3 && println("energy = $(etol/Ni/Nj*2)")
    return etol/Ni/Nj*2
end