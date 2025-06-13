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

function T_parity_conserving(T::AbstractArray)
    s = size(T)
	p = zeros(s)
	for index in CartesianIndices(T)
		i = index.I
        if sum(indextoqn.(i)) % 2 == 0 
            p[index] = 1
        end
        # i = index.I .- 1
        # sum(map((i,s)->sum(bitarray(i, ceil(Int, log2(s)))), i,s)) % 2 == 0 && (p[index] = 1)
	end
	p = _arraytype(T)(p)

	return p .* T
end

ChainRulesCore.rrule(::typeof(T_parity_conserving),T::AbstractArray) = T_parity_conserving(T), dT -> (NoTangent(), T_parity_conserving(dT))

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

function build_A(A, params::iPEPSOptimize{:fermion, :square})
    D, N = size(A)[[1,6]]
    SDD = _arraytype(A)(swapgate(D, D))
    return StructArray([T_parity_conserving(A[:,:,:,:,:,i]) for i in 1:N], params.pattern),
           StructArray([fdag(T_parity_conserving(A[:,:,:,:,:,i]), SDD) for i in 1:N], params.pattern)
end

using TeneT: mcform
function build_M(A, params::iPEPSOptimize{:fermion, :square})
    N = length(A[1])
    D = size(A[1][1], 1)
    SDD = _arraytype(A[1][1])(swapgate(D, D))
    params.ifflatten == true || throw(Base.error("ifflatten must be true for fermion currently"))
    M = StructArray([ein"((abcde,fgchi),lfbm),dkji-> glhjkema"(A[1][i],A[2][i],SDD,SDD) for i in 1:N], params.pattern)
    # Mr = reshape(M[1], D^2,D^2,D^2,D^2)
    # @show norm(M[1] - ein"(glhjkema,glbc),kedf->bcmafdhj"(conj(M[1]),SDD,SDD))
    # M = [reshape(mcform(reshape(M,D^2,D^2,D^2,D^2))[3],D,D,D,D,D,D,D,D) for M in M]
    # M = [M/norm(M) for M in M]
    return M
end
