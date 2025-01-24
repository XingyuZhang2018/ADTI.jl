using ADIT
using Test

@testset "ADIT.jl" begin
    @testset "QR factorization with $sitetype $atype{$dtype}" for atype in [Array], dtype in [ComplexF64], sitetype in [electronPn(),electronZ2(),tJZ2()]
        Random.seed!(100)
        M = randU1(sitetype, atype, dtype, 5,3,5; dir = [-1,1,1])
        Mtensor = asArray(sitetype, M)
        function foo(M)
            M = reshape(M, 15, 5)
            Q, R = qrpos(M)
            return norm(Q) + norm(R)
        end
        @test asArray(sitetype, Zygote.gradient(foo, M)[1]) ≈ num_grad(foo, Mtensor)  atol = 1e-8
    end

    @testset "LQ factorization with $sitetype $atype{$dtype}" for atype in [Array], dtype in [ComplexF64], sitetype in [electronPn(),electronZ2(),tJZ2()]
        Random.seed!(100)
        M = randU1(sitetype, atype, dtype, 5,3,5; dir = [-1,1,1])
        Mtensor = asArray(sitetype, M)
        function foo(M)
            M = reshape(M, 5, 15)
            L, Q = lqpos(M)
            return norm(Q) + norm(L)
        end
        @test asArray(sitetype, Zygote.gradient(foo, M)[1]) ≈ num_grad(foo, Mtensor)  atol = 1e-8
    end
end
