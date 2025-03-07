using LinearAlgebra
using NMRQuantumSim

# Example usage
function run_example()

    # Create a base system with 4 spins
    N_spins = 2

    # Create a base system with 4 spins
    base_system = NMRSystem(N_spins)
    set_chemical_shift!(base_system, 1, 1.5) # Fixed shift
    set_chemical_shift!(base_system, 2, -0.5) # Fixed shift

    # Define parameters to encode (J12 and h4)
    params_to_encode = ParameterSubset(
        [(1, 2)],  # Just encode J12
        []        # Just encode h4
    )

    # Create the encoding with 3 bits per parameter
    encoding = QuantumParameterEncoding(
        N_spins, 
        params_to_encode,
        j_bits_per_param=5, 
        h_bits_per_param=0,
        j_range=(0, 1.2),
        h_range=(-1, 1)
    )
    
    construct_hamiltonian_from_nmr(base_system)

    # Simulation timesteps
    ω_max = 10   # max fequency in spectrum
    γ = 0.1 # linewidth
    Δ = π / ω_max # Bandwidth relation (Nyquist)
    l = ceil(log2(2π / (Δ * γ))) # resolves linewidth
    dim = 2^N_spins
    params = NMRParameters(N_spins, l, Δ, dim, dim)

    # Compute the 4D QFT array
    A, C = compute_qft_array(base_system, encoding, params)
    println("Size of A_mnkθ: ", size(A))
    println("Size of C_mnθ: ", size(C))


    # Visualize how spectra change across parameter space
    spectral_map = visualise_parameter_space(base_system, encoding, params)
    println("Size of spectral map: ", size(spectral_map))

    fig = plot_spectral_map(spectral_map, params, "example_simulate_phase_estimation")

    return fig
    
end


@time run_example()