using LinearAlgebra

"""
    compute_A_mnk(m, n, k, eigenvalues, parameters)

Compute the A_mnk coefficient for phase estimation.
"""
function compute_A_mnk(m, n, k, eigenvalues, parameters)
    # Extract parameters
    Δ = parameters.Δ
    l = parameters.l
    
    # Energy difference - transition frequency
    ω_mn = eigenvalues[m] - eigenvalues[n]
    
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
    
    return A_mnk
end

"""
    compute_A_tensor_for_parameter(θ, parameters)

Compute the full A tensor for a given parameter set.
"""
function compute_A_tensor_for_parameter(θ, parameters)
    # Compute eigensystem for this parameter set
    eigenvalues, _ = compute_eigensystem(θ, parameters)
    
    # Extract dimensions
    num_states = length(eigenvalues)
    num_k = 2^parameters.l
    
    # Initialize result tensor
    A = Array{ComplexF64}(undef, num_states, num_states, num_k)
    
    # Compute all coefficients
    for m in 1:num_states
        for n in 1:num_states
            for k in 0:(num_k-1)
                A[m, n, k+1] = compute_A_mnk(m, n, k, eigenvalues, parameters)
            end
        end
    end
    
    return A
end

"""
    parallel_compute_A_tensors(θ_list, parameters)

Compute A tensors for multiple parameter sets.
"""
function parallel_compute_A_tensors(θ_list, parameters)
    num_θ = length(θ_list)
    
    # First compute one tensor to get dimensions
    test_tensor = compute_A_tensor_for_parameter(θ_list[1], parameters)
    num_states = size(test_tensor, 1)
    num_k = size(test_tensor, 3)
    
    # Initialize 4D result array
    A_tensors = Array{ComplexF64}(undef, num_states, num_states, num_k, num_θ)
    
    # Store the first result
    A_tensors[:, :, :, 1] = test_tensor
    
    # Compute for remaining parameter sets
    # Note: This implementation is sequential, but could be made parallel
    for i in 2:num_θ
        A_tensors[:, :, :, i] = compute_A_tensor_for_parameter(θ_list[i], parameters)
    end
    
    return A_tensors
end