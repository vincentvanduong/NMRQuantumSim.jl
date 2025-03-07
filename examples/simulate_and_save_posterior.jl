using NMRQuantumSim

# Example experiment function
function run_posterior(parameters)
    
    # Create a base system with 2 spins
    N_spins = parameters["N_spins"]
    chemical_shift1 = parameters["chemical_shift1"]
    chemical_shift2 = parameters["chemical_shift2"]

    # Create a base system with 4 spins
    base_system = NMRSystem(N_spins)
    set_chemical_shift!(base_system, 1, chemical_shift1) # Fixed shift
    set_chemical_shift!(base_system, 2, chemical_shift2) # Fixed shift

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
    n_configs = 2^total_bits(encoding)
    
    construct_hamiltonian_from_nmr(base_system)

    # Simulation timesteps
    ω_max = 10   # max fequency in spectrum
    γ = 0.05 # linewidth
    Δ = π / ω_max # Bandwidth relation (Nyquist)
    l = ceil(Int64, log2(2π / (Δ * γ))) # resolves linewidth
    dim = 2^N_spins
    n_freq = 2^l
    params = NMRParameters(N_spins, l, Δ, dim, dim)

    # Compute the 4D QFT array
    A, C = compute_qft_array(base_system, encoding, params)
    println("Size of A_mnkθ: ", size(A))
    println("Size of C_mnθ: ", size(C))

    # Make articially generate a spectrum by choosing a value of θ
    θ_ref = round(Int, n_configs / 2)
    _, A_ref = compute_spectrum(A, C, θ_ref)

    A_ref += 0.1 * randn(size(A_ref))
    
    # Calculate posterior distribution
    posterior = zeros(n_configs)

    for θ in 0:(n_configs-1)
        posterior_ = 0
        for m in 1:dim
            for n in 1:dim
                transition_amplitude = abs(C[m,n,θ+1])^2
                partial_measurement = 0
                for k in 1:n_freq-1
                    partial_measurement += A_ref[k] * A[m,n,k,θ+1]
                end
                posterior_ += transition_amplitude * abs(partial_measurement)^2
            end
        end
        posterior[θ+1] = posterior_
    end

    # In a real scenario, this would be your quantum Bayesian inference
    #estimated_params = Dict(
    #    "J_coupling" => parameters["J_coupling"] + 0.1 * randn(),
    #    "chemical_shift" => parameters["chemical_shift"] + 0.05 * randn(),
    #    "linewidth" => parameters["linewidth"] + 0.02 * randn()
    #)
    
    # Generate spectrum with estimated parameters
    #_, estimated_spectrum = generate_nmr_spectrum(estimated_params)
    
    # Calculate quality metrics
    #metrics = calculate_spectrum_quality_metrics(noisy_spectrum, estimated_spectrum)
    
    # Return results
    return Dict(
        "parameters" => collect(posterior),
        "posterior" => collect(range(0, n_configs-1, length=n_configs))
    #    "experimental_spectrum" => noisy_spectrum,
    #    "estimated_spectrum" => estimated_spectrum,
    #    "estimated_parameters" => estimated_params,
    #    "quality_metrics" => metrics
    )
end

# Run and save an experiment
parameters = Dict(
    "N_spins" => 2,             # Integer
    "J_coupling" => 3.5,        # Hz
    "chemical_shift1" => -0.5,    # ppm
    "chemical_shift2" => 1.5,    # ppm
    "linewidth" => 0.2,         # Hz
    "num_runs" => 100
)

metadata = Dict(
    "description" => "Two-spin NMR posterior example",
    "author" => "Vincent Van Duong",
    "notes" => "Simple demonstration of posterior"
)

# Run the experiment and save results
filename, results = run_and_save_experiment(
    "two_spin_simulation", 
    run_posterior, 
    parameters, 
    metadata=metadata
)
