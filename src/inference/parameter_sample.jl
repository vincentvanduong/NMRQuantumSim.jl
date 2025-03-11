# File: inference/parameter_sample.jl

"""
    ParameterSample

Structure to hold a parameter sample from the posterior distribution.

# Fields
- `system::NMRSystem`: NMR system with the sampled parameters
- `probability::Float64`: Posterior probability of this parameter configuration
- `binary_value::Int`: Binary representation of the parameters
"""
struct ParameterSample
    system::NMRSystem
    probability::Float64
    binary_value::Int
end

"""
    binary_list_to_parameter_samples(binary_values::Vector{Int}, probabilities::Vector{Float64}, 
                                    base_system::NMRSystem, encoding::QuantumParameterEncoding)

Convert a list of binary parameter configurations and probabilities to ParameterSample objects.
"""
function binary_list_to_parameter_samples(binary_values::Vector{Int}, probabilities::Vector{Float64}, 
                                        base_system::NMRSystem, encoding::QuantumParameterEncoding)
    # Check that binary_values and probabilities have the same length
    if length(binary_values) != length(probabilities)
        error("binary_values and probabilities must have the same length")
    end
    
    # Initialize result array
    samples = Vector{ParameterSample}(undef, length(binary_values))
    
    # Process each binary value
    for i in 1:length(binary_values)
        binary = binary_values[i]
        prob = probabilities[i]
        
        # Create a system with these parameters
        system = binary_to_system(binary, base_system, encoding)
        
        # Create and store ParameterSample
        samples[i] = ParameterSample(system, prob, binary)
    end
    
    return samples
end

"""
    find_map_estimate(samples::Vector{ParameterSample})

Find the maximum a posteriori parameter estimate from a collection of samples.
"""
function find_map_estimate(samples::Vector{ParameterSample})
    # Find the sample with maximum probability
    max_idx = argmax([s.probability for s in samples])
    return samples[max_idx]
end