using LinearAlgebra

"""
    create_spectral_function(H::Matrix{<:Number}, O::Matrix{<:Number}, γ::Real)

Create a callable function that efficiently computes the spectral function A(ω)
for a given Hamiltonian H, observable O, and linewidth γ.

# Arguments
- `H::Matrix{<:Number}`: Hamiltonian matrix
- `O::Matrix{<:Number}`: Observable matrix (typically total z-magnetization)
- `γ::Real`: Linewidth parameter

# Returns
- `A::Function`: Function that takes frequency ω and returns spectral intensity A(ω)
"""
function create_spectral_function(H::Matrix{<:Number}, O::Matrix{<:Number}, γ::Real)
    # Diagonalize the Hamiltonian (only once)
    eigenvalues, eigenvectors = eigen(Hermitian(H))
    
    # Calculate the matrix elements of the observable in the energy eigenbasis
    # |⟨Em|O|En⟩|²
    transition_amplitudes = Dict{Tuple{Int,Int}, Float64}()
    transition_frequencies = Dict{Tuple{Int,Int}, Float64}()
    
    # Dimension of the Hilbert space
    n_states = size(H, 1)
    
    # Compute all non-zero transition amplitudes and frequencies
    for m in 1:n_states
        for n in 1:n_states
            # Compute matrix element
            amplitude = abs(eigenvectors[:, m]' * O * eigenvectors[:, n])^2
            
            # Only store if transition is allowed (matrix element is non-zero)
            if amplitude > 1e-10
                transition_amplitudes[(m, n)] = amplitude
                transition_frequencies[(m, n)] = eigenvalues[m] - eigenvalues[n]
            end
        end
    end
    
    # Create and return a function that evaluates the spectral function at any ω
    function spectral_function(ω::Real)
        result = 0.0
        
        # Sum over all allowed transitions
        for ((m, n), amplitude) in transition_amplitudes
            ω_mn = transition_frequencies[(m, n)]
            
            # Lorentzian lineshape
            result += amplitude * γ / ((ω - ω_mn)^2 + γ^2)
        end
        
        return result
    end
    
    return spectral_function
end

"""
    create_spectral_function(system::NMRSystem, γ::Real; observable::Symbol=:z_magnetisation)

Create a spectral function for an NMRSystem object with given linewidth.

# Arguments
- `system::NMRSystem`: The NMR system to analyze
- `γ::Real`: Linewidth parameter
- `observable::Symbol=:z_magnetisation`: Type of observable to use

# Returns
- `A::Function`: Function that takes frequency ω and returns spectral intensity A(ω)
"""
function create_spectral_function(system::NMRSystem, γ::Real; observable::Symbol=:z_magnetisation)
    # Construct the Hamiltonian from the system
    H = Matrix(construct_hamiltonian_from_system(system))
    
    # Create the observable matrix
    N = system.N
    Sx, Sy, Sz = create_spin_operators(N)
    
    if observable == :z_magnetisation
        O = zeros(ComplexF64, 2^N, 2^N)
        # Total z-magnetisation S^z_tot = ∑_i S^z_i
        for i in 1:N
            O += Sz[i]
        end
    else
        error("Observable $observable not implemented")
    end
    
    return create_spectral_function(H, O, γ)
end

"""
    compute_spectrum_array(spectral_fn::Function, ω_range::AbstractRange)

Compute spectral function values over a range of frequencies.

# Arguments
- `spectral_fn::Function`: The spectral function to evaluate
- `ω_range::AbstractRange`: Range of frequencies to evaluate

# Returns
- `spectrum::Vector{Float64}`: Array of spectral intensities
"""
function compute_spectrum_array(spectral_fn::Function, ω_range::AbstractRange)
    return [spectral_fn(ω) for ω in ω_range]
end