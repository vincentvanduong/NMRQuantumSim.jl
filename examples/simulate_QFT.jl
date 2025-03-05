using NMRQuantumSim

# Example usage
function run_example()

    # Create a base system with 4 spins
    N_spins = 4
    base_system = NMRSystem(N_spins)

    # Define parameters to encode (J12 and h4)
    params_to_encode = ParameterSubset(
        [(1, 2)],  # Just encode J12
        [4]        # Just encode h4
    )

    # Create the encoding with 3 bits per parameter
    encoding = QuantumParameterEncoding(
        N_spins, 
        params_to_encode,
        j_bits_per_param=3, 
        h_bits_per_param=3,
        j_range=(-1.0, 1.0),
        h_range=(-1.0, 1.0)
    )

    # Simulation timesteps
    params = create_default_parameters(N_spins)

    # Compute the 4D QFT array
    A, C = compute_qft_array(base_system, encoding, params)
    println("Size of A_mnkθ: ", size(A))
    println("Size of C_mnθ: ", size(C))


    # Compute spectrum for parameter configuration θ=5
    frequencies, spectrum = compute_spectrum(A, C, 5)

    # Visualize how spectra change across parameter space
    spectral_map = visualize_parameter_space(base_system, encoding)
    println("Size of spectral map: ", size(spectral_map))

    
end


run_example()