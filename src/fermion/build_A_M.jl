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

function swapgate(sitetype, atype, dtype, d::Int, D::Int)
	S = ein"ij,kl->ikjl"(Matrix{dtype}(I,d,d),Matrix{dtype}(I,D,D))
	for j = 1:D, i = 1:d
        abs(indextoqn(sitetype, i)) % 2 != 0 && abs(indextoqn(sitetype, j)) % 2 != 0 && (S[i,j,i,j] = -1)
	end
	return atype(S)
end
@non_differentiable swapgate(kwarg...)

function U1swapgate(atype, dtype, d::Int, D::Int; indqn::Vector{Vector{Int}}, indims::Vector{Vector{Int}}, ifZ2::Bool)
    (d, D, d, D) != Tuple(map(sum, indims)) && throw(Base.error("U1swapgate indims is not valid"))
    dir = [-1, -1, 1, 1]
    qn = Vector{Vector{Int}}()
    tensor = Vector{atype{dtype}}()
    dims = Vector{Vector{Int}}()
    @inbounds for i in CartesianIndices(Tuple(length.(indqn)))
        qni = [indqn[j][i.I[j]] for j in 1:4]
        qnsum = ifZ2 ? sum(qni) % 2 : sum(qni .* dir)
        if qnsum == 0 && qni[1] == qni[3] && qni[2] == qni[4]
            bulkdims = [indims[j][i.I[j]] for j in 1:4]
            push!(qn, qni)
            push!(dims, bulkdims)
            tensori = atype(ein"ij,kl->ikjl"(Matrix{dtype}(I, bulkdims[1], bulkdims[1]), Matrix{dtype}(I,  bulkdims[2], bulkdims[2])))
            isodd(qni[1]) && isodd(qni[2]) && (tensori .= -tensori)
            push!(tensor, tensori)
        end
    end
    p = sortperm(qn)
    tensor = vcat(map(vec, tensor[p])...)
    return U1Array(qn[p], dir, tensor, (d, D, d, D), dims[p], 1, ifZ2)
end
@non_differentiable U1swapgate(kwarg...)

function fdag(ipeps, SDD) 
    return ein"(ulfdr,luij),pqrd->jifqp"(conj(ipeps), SDD, SDD)
end

function bulid_A(A, params::iPEPSOptimize{:fermion, :square})
    D, Ni, Nj = size(A)[[1,6,7]]
    SDD = swapgate(params.sitetype, _arraytype(A), ComplexF64, D, D)
    return [T_parity_conserving(A[:,:,:,:,:,i,j]) for i in 1:Ni, j in 1:Nj], [fdag(T_parity_conserving(A[:,:,:,:,:,i,j]), SDD) for i in 1:Ni, j in 1:Nj]
end

function bulid_M(A, params::iPEPSOptimize{:fermion, :square})
    Ni, Nj = size(A[1])
    D = size(A[1][1], 1)
    SDD = swapgate(params.sitetype, _arraytype(A[1][1]), ComplexF64, D, D)
    params.ifflatten == true || throw(Base.error("ifflatten must be true for fermion currently"))
    M = [ein"((abcde,fgchi),lfbm),dkji-> glhjkema"(A[1][i,j],A[2][i,j],SDD,SDD) for i in 1:Ni, j in 1:Nj]
    return M
end

function bulid_A(A::AbstractArray{ComplexF64, 3}, params::iPEPSOptimize{:fermion, :square})
    qnD, _, dimsD, _ = params.boundary_alg.U1info
    D = sum(dimsD)
    atype = _arraytype(A)
    Ni, Nj = size(A)[[2,3]]
    d = 4
    sitetype = params.sitetype
    info = Zygote.@ignore zerosinitial(Val(:U1), Array, ComplexF64, D,D,d,D,D; 
			dir = [-1,-1,1,1,1], 
			indqn = [qnD, qnD, getqrange(sitetype, d)..., qnD, qnD], 
			indims = [dimsD, dimsD, getblockdims(sitetype, d)..., dimsD, dimsD], 
			f=[0],
			ifZ2=sitetype.ifZ2
    )
    SDD = U1swapgate(atype, ComplexF64, D, D; 
                     indqn = [qnD for _ in 1:4], 
                     indims = [dimsD for _ in 1:4],
                     ifZ2=sitetype.ifZ2
    )
    A = [U1Array(info.qn, info.dir, atype(A[:, i, j]), info.size, info.dims, 1, sitetype.ifZ2) for i = 1:Ni, j in 1:Nj]
    Adag = [fdag(A, SDD) for A in A]
    return A, Adag
end

function bulid_M(A::Tuple{Matrix{<:U1Array}, Matrix{<:U1Array}}, params::iPEPSOptimize{:fermion, :square})
    Ni, Nj = size(A[1])
    atype = _arraytype(A[1][1].tensor)
    qnD, _, dimsD, _ = params.boundary_alg.U1info
    D = sum(dimsD)
    SDD = U1swapgate(atype, ComplexF64, D, D; 
                     indqn = [qnD for _ in 1:4], 
                     indims = [dimsD for _ in 1:4],
                     ifZ2=params.sitetype.ifZ2
    )
    params.ifflatten == true || throw(Base.error("ifflatten must be true for fermion currently"))
    M = [ein"((abcde,fgchi),lfbm),dkji-> glhjkema"(A[1][i,j],A[2][i,j],SDD,SDD) for i in 1:Ni, j in 1:Nj]
    return M
end