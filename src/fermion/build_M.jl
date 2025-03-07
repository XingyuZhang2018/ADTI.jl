function swapgate(V1::VectorSpace, V2::VectorSpace)
    I1 = isomorphism(V1, V1)
    I2 = isomorphism(V2, V2)
    S = I1 ⊗ I2
    S[(Irrep[ℤ₂](1), Irrep[ℤ₂](1), Irrep[ℤ₂](1), Irrep[ℤ₂](1))] .*= -1 # current only for Z2 symmetry 
    return S
end

@non_differentiable swapgate(kwarg...)

function fdag(ipeps::AbstractTensorMap, SDD::AbstractTensorMap) 
    @tensor ipepsdag[9 8 3 6 7] := conj(ipeps[1 2 3 4 5]) * SDD[9 8;1 2] * SDD[5 4;7 6]
    return ipepsdag
end

function build_M(A::StructArray, params::iPEPSOptimize{:fermion, :square})
    D = space(A[1], 4)
    SDD = Zygote.@ignore _arraytype(A[1])(swapgate(D, D))
    params.ifflatten == true || throw(Base.error("ifflatten must be true for fermion currently"))
    M = Zygote.Buffer(rand(ComplexF64, [D'*D ← D'*D for i in 1:length(A.data)], A.pattern))
    DD = Zygote.@ignore fuse(D', D)
    It = Zygote.@ignore isomorphism(DD, D'*D)
    for p in 1:length(A.data)
        i, j = Tuple(findfirst(==(p), A.pattern))
        Adag = fdag(A[i,j], SDD)
        @tensoropt M[i,j][-1 -2; -3 -4] := A[i,j][1 2 3 4 5] * Adag[6 7 3 8 9] * SDD[2 13; 12 6] * SDD[9 10; 11 4] * 
        conj(It[-1; 7 12]) * It[-2; 8 10] * It[-3; 11 5] * conj(It[-4; 13 1])
    end
    return copy(M)
end
