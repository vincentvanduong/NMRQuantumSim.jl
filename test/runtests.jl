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
    
    @testset "Simulation" begin
        # Test parameter creation
        params = create_default_parameters(2)
        @test params.num_spins == 2
        
        # Test Hamiltonian construction
        θ = [1.0, 2.0, 0.5]  # h1, h2, J12
        vals, vecs = compute_eigensystem(θ, params)
        @test length(vals) == 2^params.num_spins
    end
end