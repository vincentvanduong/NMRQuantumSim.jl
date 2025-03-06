using NMRQuantumSim
using LinearAlgebra

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

    fig = plot_spectral_map(spectral_map, params, "two_spin_nmr_spectral_function_with_qubits")

    return fig
    
end


"""
    compute_qft_array(base_system::NMRSystem, 
                          encoding::QuantumParameterEncoding;
                          parameters:: NMRParameters)

Compute the 4-dimensional array A_mnkθ and 3-dimensional array C_mnθ.
C_mnθ represents the scattering matrix elements of the observable with
all parameter configurations. We also compute the associated 4-dimensional array
representing the quantum Fourier transform for all parameter configurations in the encoding.

Parameters:
- base_system: The base NMR system with fixed parameters
- encoding: The parameter encoding defining which parameters vary
- parameters: The number of phase estimation qubits and the time delta.

Returns:
- A 4D array with dimensions [n_states, n_states, n_freq, n_configurations]
- A 3D array with dimensions [n_states, n_states, n_configurations]
  where n_states = 2^N for an N-spin system
"""


function compute_qft_array(base_system::NMRSystem,
    encoding::QuantumParameterEncoding,
    parameters::NMRParameters;
    observable::Symbol=:z_magnetisation)
    # Extract parameters
    Δ = parameters.Δ
    l = parameters.l
    
    N = base_system.N
    n_states = 2^N  # Hilbert space dimension
    n_configs = 2^total_bits(encoding)  # Number of parameter configurations
    n_freq = 2^l   # Number of frequencies sampled

    # Initialise 3D array
    # C[m, n, θ] represents the scattering matrix element of O from m→n for parameter θ
    C = zeros(ComplexF64, n_states, n_states, n_configs)

    Sx, Sy, Sz = create_spin_operators(N)

    if observable == :z_magnetisation
        O = zeros(ComplexF64, n_states, n_states)
        # Total z-magnetisation S^z_tot = ∑_i S^z_i
        for i in 1:N
            O += Sz[i]
        end
    else
        error("Observable $observable not implemented")
    end
    # Normalise by Tr[O²]
    normalisation = sqrt(tr(O^2))


    # Initialise the 4D array
    # A[m, n, k, θ] represents amplitude for transition m→n at frequency k for parameter θ
    A = zeros(ComplexF64, n_states, n_states, n_freq, n_configs)

    # For each parameter configuration
    for θ in 0:(n_configs-1)
        # Generate the system with these parameters
        system = binary_to_system(θ, base_system, encoding)

        # Diagonalize the Hamiltonian
        H_matrix = construct_hamiltonian_from_nmr(system)
        H_dense = Matrix(H_matrix)  # Convert sparse to dense
        eigen_decomp = eigen(Hermitian(H_dense))
        eigenvalues, eigenvectors = eigen_decomp.values, eigen_decomp.vectors

        # Compute A_mnkθ and C_mnθ for this parameter configuration
        for m in 1:n_states
            ψ_m = eigenvectors[:, m]
            for n in 1:n_states
                
                # Compute matrix element
                ψ_n = eigenvectors[:, n]
                C[m, n, θ+1] = ψ_m' * O * ψ_n

                # Transition frequency ω_mn = ε_m - ε_n
                ω_mn = eigenvalues[m] - eigenvalues[n]

                # Compute QFT coefficients for all frequency bins
                for k in 0:(n_freq-1)
                    # Frequency from phase estimation
                    ω_k = 2π * k / (Δ * 2^l)
                    
                    # Compute A_mnk coefficient using closed-form equation instead of summation
                    num = sin(Δ * 2^l * (ω_k - ω_mn) / 2)
                    denom = sin(Δ * (ω_k - ω_mn) / 2)

                    # Handle the case where denom is very close to zero
                    if abs(denom) < 1e-10
                        A_mnk = 1.0  # Use the limit value when sin(x)/x approaches 1 as x approaches 0
                    else
                        prefactor = exp(im * (2^l - 1) * Δ * (ω_k - ω_mn) / 2)
                        A_mnk = 2.0^(-l) * prefactor * num / denom
                    end
                    A[m, n, k+1, θ+1] = A_mnk
                end
            end
        end
    
    end

    C ./= normalisation

    return (A, C)
end


"""
    compute_spectrum(A::Array{ComplexF64, 4}, C::Matrix{ComplexF64}, θ::Integer; 
                    γ::Float64=0.5, normalise::Bool=true)

Compute the NMR spectrum for a specific parameter configuration θ.

Parameters:
- A: The 4D array from compute_qft_array
- C: The 3D transition amplitude matrix
- θ: The parameter configuration index (0-based)
- γ: Line broadening factor
- normalise: Whether to normalise the spectrum

Returns:
- frequencies: Array of frequency points
- spectrum: Array of spectral intensities
"""
function compute_spectrum(A::Array{ComplexF64, 4}, C::Array{ComplexF64, 3}, θ::Integer; 
                         normalise::Bool=true)
    
    n_states, _, n_freq, _ = size(A)
    
    # Extract the A[:,:,:,θ+1] and C[:,:, θ+1] slice
    A_θ = view(A, :, :, :, θ+1)
    C_θ = view(C, :, :, θ+1)
    
    # Compute spectrum
    spectrum = zeros(Float64, n_freq)
    for k in 1:n_freq
        for m in 1:n_states
            for n in 1:n_states
                # Add contribution to spectrum
                spectrum[k] += abs(C_θ[m,n])^2 * (abs(A_θ[m,n,k]))^2
            end
        end
    end
    
    # normalise if requested
    if normalise && maximum(spectrum) > 0
        spectrum ./= maximum(spectrum)
    end
    
    # Generate frequency points (Normalised to [0,1])
    frequencies = collect(range(0, 1, length=n_freq))
    
    return frequencies, spectrum
end


function visualise_parameter_space(base_system::NMRSystem,
    encoding::QuantumParameterEncoding,
    parameters::NMRParameters)

    resolution = 2^(parameters.l)

    
    # Compute the QFT array
    A, C = compute_qft_array(base_system, encoding, parameters)
    
    # Number of parameter configurations
    n_configs = 2^total_bits(encoding)
    
    # Initialize spectral map
    spectral_map = zeros(n_configs, resolution)
    
    # Compute spectrum for each parameter configuration
    for θ in 0:(n_configs-1)
        _, spectrum = compute_spectrum(A, C, θ)
        spectral_map[θ+1, :] = spectrum
    end
    
    return spectral_map
end





@time run_example()