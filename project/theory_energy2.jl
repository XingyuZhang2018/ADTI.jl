using QuadGK

# 被积函数定义
function integrand(kx, ky)
    d = √2 * (cos(kx) - cos(ky))
    c1 = exp(im*π/4) * (1 + exp(-im*(kx + ky)))
    c2 = exp(-im*π/4) * (exp(-im*kx) + exp(-im*ky))
    c = c1 + c2
    return real(-√(abs2(d) + abs2(c)))  # 取实部保证数值稳定性
end

# 双重数值积分
result, err = quadgk(kx -> quadgk(ky -> integrand(kx, ky), 0, 2π)[1], 0, 2π, rtol=1e-9)
average_energy = result / (4π^2)

println("基态平均能量: ", average_energy)
println("积分误差估计: ", err)