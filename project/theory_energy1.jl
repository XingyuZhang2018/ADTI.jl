using LinearAlgebra

# 定义常数
t1 = exp(im * π / 4)
t1c = conj(t1)
t2 = 1 / √2

# 定义哈密顿量函数
function H(kx::Float64, ky::Float64)
    H_matrix = zeros(ComplexF64, 2, 2)
    H_matrix[1, 2] = t1 + t1c * exp(-im * kx) + t1c * exp(-im * ky) + t1 * exp(-im * (kx + ky))
    H_matrix[2, 1] = conj(H_matrix[1, 2])
    diag_term = 2 * t2 * (cos(kx) - cos(ky))
    H_matrix[1, 1] = diag_term
    H_matrix[2, 2] = -diag_term
    return H_matrix
end

# 参数设置
Lx = 1000
Ly = 1000
kxlist = 2π / Lx .* (0:Lx-1)
kylist = 2π / Ly .* (0:Ly-1)

# 初始化能量数组
Elist = zeros(2, Lx, Ly)

# 计算能带
for i in 1:Lx
    for j in 1:Ly
        kx = kxlist[i]
        ky = kylist[j]
        H_matrix = H(kx, ky)
        e = eigvals(Hermitian(H_matrix))  # 计算厄米矩阵的本征值（升序排列）
        Elist[1, i, j] = e[end]  # 最高能带
        Elist[2, i, j] = e[1]    # 最低能带
    end
end

# 计算最低能带的平均值
result = sum(Elist[2, :, :]) / (Lx * Ly)
println("Average energy of the lower band: ", result)