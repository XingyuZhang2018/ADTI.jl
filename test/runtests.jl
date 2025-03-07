using ADTI
using ADTI: hamiltonian, const_Sx, const_Sy, const_Sz
using ADTI: swapgate
using ADTI: build_M
using Random
using Test
using TeneT
using TensorKit
using OptimKit
using OMEinsum
@testset "ADTI.jl" begin


    @testset "test_hamiltonian_models" begin
        @testset "H trunc" begin
            Sx = const_Sx(1.0)
            Sy = const_Sy(1.0)
            Sz = const_Sz(1.0)
            H = ein"ij,kl->ijkl"(Sx, Sx) + ein"ij,kl->ijkl"(Sy, Sy) + ein"ij,kl->ijkl"(Sz, Sz)
            U, S, V = svd(reshape(H, 9,9))
            @test sum(S.>1e-10) == 3
        end
        @testset "Hubbard" begin
            model = Hubbard(1.0, 12.0, 6.0)
            h1, h2 = hamiltonian(model)
            H = h1 * h2
            H_array = reshape(permutedims(convert(Array, H), (1,3,2,4)), 4^2, 4^2)
            @test H_array ≈ H_array'
            U, S, Vd = tsvd(H, trunc = truncbelow(1e-10))
            @test dim(S) == 4*4+2*2
            @show S
            A = U * sqrt(S)
            B = sqrt(S) * Vd
            @test H ≈ A * B 
        end
        @testset "Topological_Insulator" begin
            model = Topological_Insulator(1.0, pi/4)
            h1, h2 = hamiltonian(model)
            H = h1 * h2
            H_array = reshape(permutedims(convert(Array, H), (1,3,2,4)), 4^2, 4^2)
            @test H_array ≈ H_array'
            U, S, Vd = tsvd(H, trunc = truncbelow(1e-10))
            @test dim(S) == 4*4
            A = U * sqrt(S)
            B = sqrt(S) * Vd
            @test H ≈ A * B 
        end
    end
    @testset "init_ipeps" begin
        system = :fermion
        lattice = :square
        boundary_alg = VUMPS()
        pattern = [1 2; 2 1]
        params = iPEPSOptimize{system, lattice}(boundary_alg=boundary_alg, verbosity=0)
        d = ℤ₂Space(0=>2, 1=>2)
        D = ℤ₂Space(0=>2, 1=>2)
        χ = ℤ₂Space(0=>2, 1=>2)
        A = init_ipeps(;atype = Array, params, d, D, χ, pattern)
        @test size(A) == (2,2)
        @test A[1,1] == A[2,2]
        @test A[1,2] == A[2,1]
        @test A[1].space.codomain == D'*D'*d*D*D
    end

    @testset "swapgate" begin
        V1 = ℤ₂Space(0=>1, 1=>1)
        V2 = ℤ₂Space(0=>1, 1=>1)
        S = swapgate(V1, V2)
        @test S.space == (V1*V2 ← V1*V2)
    end

    @testset "build_M" begin
        system = :fermion
        lattice = :square
        boundary_alg = VUMPS()
        pattern = [1 2; 2 1]
        params = iPEPSOptimize{system, lattice}(boundary_alg=boundary_alg, verbosity=0)
        d = ℤ₂Space(0=>1, 1=>1)
        D = ℤ₂Space(0=>2, 1=>2)
        χ = ℤ₂Space(0=>3, 1=>3)
        @test typeof(d) <: VectorSpace
        A = init_ipeps(;atype = Array, params, No=0, d, D, χ, pattern)

        M = build_M(A, params)
        @test size(M) == (2,2)
        @test M[1,1] == M[2,2]
        @test M[1,2] == M[2,1]
        @test M[1].space == (fuse(D'*D)'*fuse(D'*D) ← fuse(D'*D)'*fuse(D'*D))
    end

    @testset "energy for pattern = $pattern" for pattern in [[1;;], [1 2; 2 1]]
        Random.seed!(100)
        system = :fermion
        lattice = :square
        boundary_alg = VUMPS(verbosity=0, ifupdown=true, maxiter=10, ifdownfromup=true)
        params = iPEPSOptimize{system, lattice}(boundary_alg=boundary_alg, verbosity=0)
        d = ℤ₂Space(0=>2, 1=>2)
        D = ℤ₂Space(0=>1, 1=>1)
        χ = ℤ₂Space(0=>10, 1=>10)
        A = init_ipeps(;atype = Array, params, No=0, d, D, χ, pattern)
        model = Free_Fermion()
        M = build_M(A, params)
        rt = VUMPSRuntime(M, χ, params.boundary_alg)
        @test imag(energy(A, model, rt, params)) < 1e-7
    end

    @testset "optimise_ipeps" begin
        Random.seed!(100)
        system = :fermion
        lattice = :square
        boundary_alg = VUMPS(verbosity=0, ifupdown=true, maxiter=10, ifdownfromup=true)
        pattern = [1;;]
        params = iPEPSOptimize{system, lattice}(boundary_alg=boundary_alg, 
                                                optimizer=LBFGS(; gradtol=1e-6, maxiter=1),
                                                verbosity=0)
        d = ℤ₂Space(0=>2, 1=>2)
        D = ℤ₂Space(0=>1, 1=>1)
        χ = ℤ₂Space(0=>10, 1=>10)
        A = init_ipeps(;atype = Array, params, No=0, d, D, χ, pattern)
        model = Free_Fermion()
        x, f, g, numfg, normgradhistory = optimise_ipeps(A, model, χ, params)
        @test f ≈ -0.4681404325943741
    end
end;
