# File: test/test_spectrum.jl

using Test
using NMRQuantumSim

@testset "Basic specrtum evaluation" begin
    system = NMRSystem(4)  # 4-spin system
    set_coupling!(system, 1, 2, 0.5)
    set_coupling!(system, 2, 3, 0.3)
    set_chemical_shift!(system, 1, 0.1)
    set_chemical_shift!(system, 2, -0.2)

    # Create a spectral function with linewidth γ = 0.01
    γ = 0.01
    spectral_fn = create_spectral_function(system, γ)

    # Evaluate at a single frequency
    ω = 0.5
    A_value = spectral_fn(ω)
    @test A_value > 0
    println("Evaluate A(ω) at ω=$ω: ", A_value)
end