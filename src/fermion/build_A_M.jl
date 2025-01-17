function indextoqn(i::Int)
    i -= 1
    i == 0 && return 0

    qni = []
    while i > 0
        remainder = i % 4
        if remainder == 2
            remainder = 1
        elseif remainder == 3
            remainder = 0
        end
        pushfirst!(qni, remainder)
        i = div(i, 4)
        # remainder = i % 2
        # pushfirst!(qni, remainder)
        # i = div(i, 2)
    end

    return sum(qni) % 2
end

function swapgate(d::Int, D::Int)
	S = ein"ij,kl->ikjl"(Matrix{ComplexF64}(I,d,d),Matrix{ComplexF64}(I,D,D))
	for j = 1:D, i = 1:d
        if abs(indextoqn(i)) % 2 != 0 && abs(indextoqn(j)) % 2 != 0
            S[i,j,i,j] = -1
        end
	end
	return S
end

@non_differentiable swapgate(d, D)

function fdag(ipeps, SDD) 
    return ein"(ulfdr,luij),pqrd->jifqp"(conj(ipeps), SDD, SDD)
end

function bulid_A(A, ::iPEPSOptimize{:fermion, :square})
    D, Ni, Nj = size(A)[[1,6,7]]
    SDD = _arraytype(A)(swapgate(D, D))
    return [A[:,:,:,:,:,i,j] for i in 1:Ni, j in 1:Nj], [fdag(A[:,:,:,:,:,i,j], SDD) for i in 1:Ni, j in 1:Nj]
end

function bulid_M(A, params::iPEPSOptimize{:fermion, :square})
    Ni, Nj = size(A[1])
    D = size(A[1][1], 1)
    SDD = _arraytype(A[1])(swapgate(D, D))
    params.ifflatten == true || throw(Base.error("ifflatten must be true"))
    M = [reshape(ein"((abcde,fgchi),lfbm),dkji-> glhjkema"(A[1][i,j],A[2][i,j],SDD,SDD), D^2,D^2,D^2,D^2) for i in 1:Ni, j in 1:Nj]
    return M
end
