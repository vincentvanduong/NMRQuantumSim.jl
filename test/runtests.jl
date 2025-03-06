using Test
using NMRQuantumSim

@testset "NMRQuantumSim.jl" begin
    @testset "NMRSystem" begin
        # Basic system creation test
        system = NMRSystem(3)
        @test system.N == 3
        @test size(system.J) == (3, 3)
        @test length(system.h) == 3
    end
    
end