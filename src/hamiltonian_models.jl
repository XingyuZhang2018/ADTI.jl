function const_Sx(S::Real)
    dims = Int(2*S + 1)
    ms = [-S+i-1 for i in 1:dims]
    Sx = zeros(ComplexF64, dims, dims)
    for j in 1:dims, i in 1:dims
        if abs(i-j) == 1
            Sx[i,j] = 1/2 * sqrt(S*(S+1)-ms[i]*ms[j]) 
        end
    end
    return Sx
end

function const_Sy(S::Real)
    dims = Int(2*S + 1)
    ms = [-S+i-1 for i in 1:dims]
    Sy = zeros(ComplexF64, dims, dims)
    for j in 1:dims, i in 1:dims
        if i-j == 1
            Sy[i,j] = -1/2/1im * sqrt(S*(S+1)-ms[i]*ms[j]) 
        elseif j-i == 1
            Sy[i,j] =  1/2/1im * sqrt(S*(S+1)-ms[i]*ms[j]) 
        end
    end
    return Sy
end

function const_Sz(S::Real)
    dims = Int(2*S + 1)
    ms = [S-i+1 for i in 1:dims]
    Sz = zeros(ComplexF64, dims, dims)
    for i in 1:dims
        Sz[i,i] = ms[i]
    end
    return Sz
end

abstract type HamiltonianModel end
@doc raw"
    Kitaev(Jx::T,Jz::T,Jy::T) where {T<:Real}

return a struct representing the Kitaev model with magnetisation fields
`Jx`, `Jy` and `Jz`..
"
struct Kitaev <: HamiltonianModel
    S::Real
    Jx::Real
    Jy::Real
    Jz::Real
end
Kitaev() = Kitaev(0.5, -1.0, -1.0, -1.0)

struct Hubbard <: HamiltonianModel
    t::Real
	U::Real
    μ::Real
end

function hamiltonian(model::Hubbard)
    t = model.t
	U = model.U
    μ = model.μ
    H = zeros(4,4,4,4)

    # t
    H[1,2,2,1], H[1,4,2,3], H[2,1,1,2], H[2,3,1,4] = -t, -t, -t, -t
    H[1,3,3,1], H[2,3,4,1], H[3,1,1,3], H[4,1,2,3] = -t, -t, -t, -t
    H[3,2,4,1], H[3,4,4,3], H[4,1,3,2], H[4,3,3,4] = t, t, t, t
    H[1,4,3,2], H[2,4,4,2], H[3,2,1,4], H[4,2,2,4] = t, t, t, t

    # U
    H[1,4,1,4], H[4,1,4,1] = U/4, U/4
    H[2,4,2,4], H[4,2,4,2] = U/4, U/4
    H[3,4,3,4], H[4,3,4,3] = U/4, U/4
    H[4,4,4,4] = U/2

    # μ
    H[1,2,1,2], H[2,1,2,1], H[1,3,1,3], H[3,1,3,1] = -μ/4, -μ/4, -μ/4, -μ/4
    H[2,3,2,3], H[3,2,3,2], H[2,2,2,2], H[3,3,3,3] = -μ/2, -μ/2, -μ/2, -μ/2
    H[1,4,1,4] += -μ/2
    H[4,1,4,1] += -μ/2
    H[2,4,2,4] += -3*μ/4
    H[4,2,4,2] += -3*μ/4
    H[3,4,3,4] += -3*μ/4
    H[4,3,4,3] += -3*μ/4
    H[4,4,4,4] += -μ

    return H
end
@non_differentiable hamiltonian(model)