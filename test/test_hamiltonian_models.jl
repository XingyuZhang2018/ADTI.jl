@testset "Topological_Insulator" begin
    model = Topological_Insulator(1.0, pi/4)
    H = hamiltonian(model)
    @test size(H) == (4,4,4,4)
    Hm = reshape(H, 16, 16)
    @test Hm ≈ Hm'
end