"""
    compute_overlap_amplitude(target_spectrum::Vector{Float64}, A::Array{ComplexF64, 4}, 
                             C::Array{ComplexF64, 3})

Compute the amplitude of overlap between the target spectrum and the model spectrum for all 
parameter configurations in the encoding.
"""
function compute_likelihood(target_spectrum::Union{Vector{ComplexF64},Vector{Float64}},
                                  A::Array{ComplexF64, 4}, 
                                  C::Array{ComplexF64, 3})
    n_states, _, n_freq, n_configs = size(A)
    
    # Pre-normalize target spectrum
    norm_target = target_spectrum / sqrt(sum(abs2, target_spectrum))
    
    # Initialize likelihood for each configuration
    likelihoods = zeros(Float64, n_configs)
    
    # Compute overlap for each parameter configuration
    for θ in 1:n_configs
        # Compute model spectrum for this configuration
        amplitude = 0
        for m in 1:n_states
            for n in 1:n_states:n_states
                overlap = 0
                for k in 1:n_freq
                    overlap += conj(norm_target[k]) * A[m,n,k,θ]
                end
                amplitude += abs2(C[m,n,θ]) * abs2(overlap)
            end
        end
        likelihoods[θ] = amplitude
    end

    binary_values = collect(0:(n_configs-1))
    
    return likelihoods, binary_values
end


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
    posterior = likelihood .* prior
    
    # Normalize posterior
    posterior ./= sum(posterior)
    
    # Convert to parameter samples
    samples = binary_list_to_parameter_samples(binary_values, posterior, base_system, encoding)
    
    return posterior, samples
end