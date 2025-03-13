# File: test/test_posterior_simple.jl

using Test
using NMRQuantumSim

@testset "Basic Posterior Estimation" begin
    # Create a base system with 2 spins
    N_spins = 2
    base_system = NMRSystem(N_spins)
    set_coupling!(base_system, 1, 2, 0.4)      # Variable coupling
    set_chemical_shift!(base_system, 1, 0.5)  # Variable shift
    set_chemical_shift!(base_system, 2, -1.5) # Fixed shift

    # Define which parameters to encode
    params_to_encode = ParameterSubset(
        [(1, 2)],  # Just encode J12
        []        
    )

    # Create the encoding with 3 bits per parameter
    encoding = QuantumParameterEncoding(
        N_spins, # Number of spins
        params_to_encode,
        j_bits_per_param=8, 
        h_bits_per_param=0,
        j_range=(0, 1.0),
        h_range=(0, 1.0)
    )

    # Calculate total qubits needed
    n_qubits = total_bits(encoding)
    @test n_qubits == 8  # Should be 10 (10 for J12)

    # Get parameter-to-qubit mapping
    param_map = parameter_index_map(encoding)
    println("J[1,2] is encoded in qubits: ", param_map[(:J, 1, 2)])

    # Let's test with a specific parameter setting
    true_system = deepcopy(base_system)
    true_parameters = system_to_binary(true_system, encoding)

    # Convert to binary and back
    binary_params = system_to_binary(true_system, encoding)
    println("Binary representation: ", bitstring(binary_params))

    # Decode and verify
    decoded_system = binary_to_system(binary_params, base_system, encoding)
    println("Original J[1,2]: ", true_system.J[1,2], ", Decoded J[1,2]: ", decoded_system.J[1,2])

    @test true_system.J[1,2] ≈ decoded_system.J[1,2] atol=0.1

    # Simulation timesteps
    ω_max = 10   # max fequency in spectrum
    γ = 0.05 # linewidth
    Δ = π / ω_max # Bandwidth relation (Nyquist)
    l = ceil(Int64, log2(2π / (Δ * γ))) # resolves linewidth
    dim = 2^N_spins
    n_freq = 2^l
    parameters = NMRParameters(N_spins, l, Δ, dim, dim)

    # Compute the target spectrum
    A, C = compute_qft_array(base_system, encoding, parameters)
    target_spectrum = compute_spectrum(A, C, true_parameters)
    #target_spectrum += 0.05 * abs.(randn(size(target_spectrum))) # Add a bit of noise

    samples = estimate_posterior(target_spectrum, base_system, encoding, parameters)
    best_config = find_map_estimate(samples)

    println(samples)
    
    # Check if recovered parameters are close to the test values
    println("Original J[1,2]: ", true_system.J[1,2], ", Best J[1,2]: ", best_config.parameter_values[:J_1_2])
    @test isapprox(true_system.J[1,2], best_config.parameter_values[:J_1_2], atol=0.1)
    

end