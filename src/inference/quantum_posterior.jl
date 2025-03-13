"""
    compute_likelihood(target_spectrum::Vector{Float64}, A::Array{ComplexF64, 4}, 
                      C::Array{ComplexF64, 3})

Compute the likelihood of the target spectrum for all parameter configurations in the encoding.
This quantifies how well each parameter configuration explains the observed spectrum.
"""
function compute_likelihood(target_spectrum::Union{Vector{ComplexF64},Vector{Float64}},
                           A::Array{ComplexF64, 4}, 
                           C::Array{ComplexF64, 3})
    n_states, _, n_freq, n_configs = size(A)
    
    # Pre-normalize target spectrum
    norm_target = target_spectrum / sqrt(sum(abs2, target_spectrum))
    
    # Pre-compute C squared magnitudes
    C_squared = abs2.(C)
    
    # Initialize likelihood for each configuration
    likelihoods = zeros(Float64, n_configs)
    
    # Compute overlap for each parameter configuration
    for θ in 1:n_configs
        # Reshape target spectrum for broadcasting
        norm_target_reshaped = reshape(conj.(norm_target), 1, 1, n_freq)
        
        # Compute overlaps for all m,n pairs at once
        overlaps = sum(norm_target_reshaped .* view(A, :, :, :, θ), dims=3)
        
        # Square the magnitudes
        overlaps_squared = abs2.(overlaps)
        
        # Multiply with C_squared and sum
        likelihoods[θ] = sum(C_squared[:,:,θ] .* overlaps_squared)
    end

    binary_values = collect(0:(n_configs-1))
    
    return likelihoods, binary_values
end


"""
    estimate_posterior(target_spectrum, base_system, encoding, parameters; prior=nothing)

Estimate the posterior distribution over model parameters given an observed target spectrum.
Returns samples from the posterior distribution.

Parameters:
- target_spectrum: The experimentally observed spectrum
- base_system: The NMR system with fixed parameters
- encoding: The parameter encoding scheme
- parameters: Control parameters for phase estimation
- prior: Optional prior distribution (uniform if not specified)

Returns:
- samples: Collection of posterior samples with parameter values and probabilities
"""
function estimate_posterior(target_spectrum::Union{Vector{ComplexF64},Vector{Float64}},
                           base_system::NMRSystem,
                           encoding::QuantumParameterEncoding,
                           parameters::NMRParameters;
                           prior::Union{Vector{Float64},Nothing}=nothing
                           )
    
    # Compute phase estimation arrays
    A, C = compute_qft_array(base_system, encoding, parameters)
    
    # Compute likelihood
    likelihood, binary_values = compute_likelihood(target_spectrum, A, C)
    
    # Get number of configurations
    n_configs = length(likelihood)
    
    # Use uniform prior if none provided
    if isnothing(prior)
        prior = ones(Float64, n_configs) ./ n_configs
    end
    
    # Apply Bayes' rule: posterior ∝ likelihood × prior
    # Use element-wise multiplication (.* operator) for vectors
    posterior = likelihood .* prior
    
    # Normalize posterior
    posterior = posterior ./ sum(posterior)  # Fixed: Changed 'probability' to 'posterior'
    
    # Convert to parameter samples
    samples = collect_posterior_samples(binary_values, posterior, base_system, encoding)
    
    return samples
end