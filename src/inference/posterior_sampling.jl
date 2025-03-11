"""
    BinaryPosteriorSample

Structure to hold a single sample from the posterior distribution in both binary and parameter space.
"""
struct BinaryPosteriorSample
    binary_value::Int
    parameter_values::Dict{Symbol,Float64}
    log_probability::Float64
end

"""
    convert_binary_sample_to_parameters(binary_sample::Integer, base_system::NMRSystem, 
                                       encoding::QuantumParameterEncoding)

Convert a binary sample to a dictionary of physical parameter values.
"""
function convert_binary_sample_to_parameters(binary_sample::Integer, base_system::NMRSystem, 
                                           encoding::QuantumParameterEncoding)
    # Create a modified system with the sampled parameters
    modified_system = binary_to_system(binary_sample, base_system, encoding)
    
    # Extract the parameters into a dictionary
    params = Dict{Symbol,Float64}()
    
    # Extract J couplings
    for (i, j) in encoding.params.j_couplings
        param_key = Symbol("J_$(i)_$(j)")
        params[param_key] = modified_system.J[i,j]
    end
    
    # Extract chemical shifts
    for i in encoding.params.h_shifts
        param_key = Symbol("h_$(i)")
        params[param_key] = modified_system.h[i]
    end
    
    return params
end

"""
    collect_posterior_samples(binary_samples::Vector{Int}, log_probs::Vector{Float64}, 
                             base_system::NMRSystem, encoding::QuantumParameterEncoding)

Process raw binary samples from quantum circuit into parameter space samples.
"""
function collect_posterior_samples(binary_samples::Vector{Int}, log_probs::Vector{Float64}, 
                                  base_system::NMRSystem, encoding::QuantumParameterEncoding)
    samples = Vector{BinaryPosteriorSample}(undef, length(binary_samples))
    
    for i in 1:length(binary_samples)
        binary_value = binary_samples[i]
        params = convert_binary_sample_to_parameters(binary_value, base_system, encoding)
        samples[i] = BinaryPosteriorSample(binary_value, params, log_probs[i])
    end
    
    return samples
end

"""
    estimate_posterior_density(samples::Vector{BinaryPosteriorSample}, 
                              param_names::Vector{Symbol}; 
                              bandwidth=:silverman)

Estimate the posterior probability density using kernel density estimation.
"""
function estimate_posterior_density(samples::Vector{BinaryPosteriorSample}, 
                                   param_names::Vector{Symbol}; 
                                   bandwidth=:silverman)
    # Extract parameter values into a matrix
    n_samples = length(samples)
    n_params = length(param_names)
    
    X = zeros(n_samples, n_params)
    for i in 1:n_samples
        for (j, param) in enumerate(param_names)
            X[i, j] = samples[i].parameter_values[param]
        end
    end
    
    # Calculate weights from log probabilities
    log_probs = [s.log_probability for s in samples]
    # Normalize log probabilities to avoid numerical issues
    log_probs .-= maximum(log_probs)
    probs = exp.(log_probs)
    weights = probs ./ sum(probs)
    
    # Construct kernel density estimate
    # Note: You would need a KDE package like KernelDensity.jl
    # kde = kde_estimate(X, weights, bandwidth)
    
    # For demonstration, we'll return the raw data
    return (X, weights, param_names)
end

"""
    marginal_posterior(density_estimate, param_index::Int; grid_size::Int=100)

Extract a marginal posterior distribution for a single parameter.
"""
function marginal_posterior(density_estimate, param_index::Int; grid_size::Int=100)
    X, weights, param_names = density_estimate
    
    # Extract the values for the specified parameter
    param_values = X[:, param_index]
    param_name = param_names[param_index]
    
    # Create a grid for evaluation
    min_val, max_val = extrema(param_values)
    grid = range(min_val, max_val, length=grid_size)
    
    # For demonstration, we'll create a simple histogram
    hist = fit(Histogram, param_values, weights=weights, nbins=grid_size)
    
    return (param_name, grid, hist.weights ./ sum(hist.weights))
end

"""
    find_map_estimate(samples::Vector{BinaryPosteriorSample})

Find the maximum a posteriori estimate from the samples.
"""
function find_map_estimate(samples::Vector{BinaryPosteriorSample})
    # Find the sample with the highest posterior probability
    max_idx = argmax([s.log_probability for s in samples])
    return samples[max_idx]
end

"""
    parameter_credible_intervals(samples::Vector{BinaryPosteriorSample}, 
                                param_names::Vector{Symbol};
                                alpha::Float64=0.05)

Calculate credible intervals for each parameter.
"""
function parameter_credible_intervals(samples::Vector{BinaryPosteriorSample}, 
                                     param_names::Vector{Symbol};
                                     alpha::Float64=0.05)
    n_samples = length(samples)
    intervals = Dict{Symbol, Tuple{Float64, Float64}}()
    
    for param in param_names
        # Extract values for this parameter
        values = [s.parameter_values[param] for s in samples]
        
        # Sort values
        sort!(values)
        
        # Calculate credible interval indices
        lower_idx = ceil(Int, alpha/2 * n_samples)
        upper_idx = floor(Int, (1-alpha/2) * n_samples)
        
        # Store interval
        intervals[param] = (values[lower_idx], values[upper_idx])
    end
    
    return intervals
end

