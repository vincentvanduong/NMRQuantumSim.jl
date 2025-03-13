using LinearAlgebra
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
        H = Matrix(construct_hamiltonian_from_system(system))
        eigenvalues, eigenvectors = compute_eigensystem(H)

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
                    #TODO Handle case where the limit value approaches zero.
                    phase_diff = (ω_k - ω_mn)
                    
                    # Compute A_mnk coefficient using closed-form equation instead of summation
                    num = sin(Δ * 2^l * (ω_k - ω_mn) / 2)
                    denom = sin(Δ * (ω_k - ω_mn) / 2)

                    # Handle the case where denom is very close to zero
                    if abs(phase_diff) < 1e-10
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
    compute_spectrum(A::Array{ComplexF64, 4}, C::Array{ComplexF64, 3}, θ::Integer; 
                     normalise::Bool=false)

Compute the NMR spectrum for a specific parameter configuration θ.

# Arguments
- `A::Array{ComplexF64, 4}`: The 4D array containing frequency components from quantum Fourier transform
- `C::Array{ComplexF64, 3}`: The 3D transition amplitude matrix
- `θ::Integer`: The parameter configuration index (0-based)
- `normalise::Bool=false`: Whether to normalize the spectrum

# Returns
- `spectrum::Vector{Float64}`: Array of spectral intensities
"""
function compute_spectrum(A::Array{ComplexF64, 4}, C::Array{ComplexF64, 3}, θ::Integer; 
                         normalise::Bool=false)
    
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
    
    return spectrum
end

"""
    visualise_parameter_space(base_system::NMRSystem,
                             encoding::QuantumParameterEncoding,
                             parameters::NMRParameters)

Generate a visualization of the NMR spectra across all parameter configurations in the encoded parameter space.

# Arguments
- `base_system::NMRSystem`: The quantum system representing the NMR setup
- `encoding::QuantumParameterEncoding`: The parameter encoding scheme
- `parameters::NMRParameters`: Parameters controlling the phase estimation and spectral resolution

# Returns
- `spectral_map::Matrix{Float64}`: A matrix where each row contains the spectrum for a specific parameter configuration
"""
function visualise_parameter_space(base_system::NMRSystem,
    encoding::QuantumParameterEncoding,
    parameters::NMRParameters)

    n_freq = 2^parameters.l

    # Compute the QFT array
    A, C = compute_qft_array(base_system, encoding, parameters)
    
    # Number of parameter configurations
    n_configs = 2^total_bits(encoding)
    
    # Initialize spectral map
    spectral_map = zeros(n_configs, n_freq)
    
    # Compute spectrum for each parameter configuration
    for θ in 0:(n_configs-1)
        # Fix: compute_spectrum returns a single value, not a tuple
        spectrum = compute_spectrum(A, C, θ)
        spectral_map[θ+1, :] = spectrum
    end
    
    return spectral_map
end