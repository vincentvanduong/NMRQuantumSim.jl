using LinearAlgebra

"""
    NMRSystem{T<:AbstractFloat}

A struct representing an NMR spin system with its Hamiltonian parameters.
"""

struct NMRSystem
    N::Int
    J::Matrix{Float64}
    h::Vector{Float64}
    labels::Vector{String}
    
    # Constructor with validation
    function NMRSystem(N::Int, J::Matrix{Float64}, h::Vector{Float64}, 
                       labels::Vector{String}=["spin$i" for i in 1:N])
        # Validate dimensions
        size(J) == (N, N) || error("J matrix must be $N × $N")
        length(h) == N || error("h vector must have length $N")
        length(labels) == N || error("labels vector must have length $N")
        
        # Ensure J is symmetric (J_{ij} = J_{ji})
        issymmetric(J) || error("J matrix must be symmetric")
        
        # Ensure diagonal of J is zero (no self-coupling)
        all(diag(J) .== 0.0) || error("Diagonal elements of J must be zero")
        
        new(N, J, h, labels)
    end
end

# Convenience constructors

# Create an NMRSystem with N spins and all parameters initialized to zero.
function NMRSystem(N::Int)
    NMRSystem(N, zeros(N, N), zeros(N))
end

# Only chemical shifts (no J coupling)
function NMRSystem(h::Vector{Float64})
    N = length(h)
    NMRSystem(N, zeros(N, N), h)
end

# Random system
function initialize_random_system(N::Int; 
                                j_range::Tuple{Float64,Float64}=(-10.0, 10.0),
                                h_range::Tuple{Float64,Float64}=(0.0, 100.0))
    # Create empty J matrix
    J = zeros(N, N)
    
    # Fill upper triangular part (excluding diagonal)
    for i in 1:N
        for j in (i+1):N
            J[i,j] = rand() * (j_range[2] - j_range[1]) + j_range[1]
        end
    end
    
    # Make it symmetric
    J = J + J'
    
    # Generate random chemical shifts
    h = rand(N) .* (h_range[2] - h_range[1]) .+ h_range[1]
    
    NMRSystem(N, J, h)
end

# Parameter manipulation functions

# Modify the coupling between sites i and j
function set_coupling!(system::NMRSystem, i::Int, j::Int, value::Float64)
    system.J[i,j] = value
    system.J[j,i] = value  # Maintain symmetry
    return system
end

# Modify the chemical shift at site i
function set_chemical_shift!(system::NMRSystem, i::Int, value::Float64)
    system.h[i] = value
    return system
end

# Extract all parameters from the system as a flat vector.
function hamiltonian_parameters(system::NMRSystem)
    # Extract upper triangular part of J (excluding diagonal)
    j_params = [system.J[i,j] for i in 1:system.N for j in (i+1):system.N]
    
    # Concatenate with chemical shifts
    return vcat(j_params, system.h)
end

# Update the system with new parameters
function update_from_parameters!(system::NMRSystem, params::Vector{Float64})
    N = system.N
    n_j_params = N * (N - 1) ÷ 2  # Number of J parameters (upper triangular)
    
    # Check if the parameter vector has the correct length
    if length(params) != n_j_params + N
        error("Parameter vector must have length $(n_j_params + N)")
    end
    
    # Update J parameters (upper triangular part)
    param_idx = 1
    for i in 1:N
        for j in (i+1):N
            system.J[i,j] = system.J[j,i] = params[param_idx]
            param_idx += 1
        end
    end
    
    # Update chemical shifts
    system.h .= params[n_j_params+1:end]
    
    return system
end