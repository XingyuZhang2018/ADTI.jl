using ADTI
using ADTI: hamiltonian
using Test

@testset "ADTI.jl" begin
    include("test_hamiltonian_models.jl")
end
