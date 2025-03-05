using StaticArrays

"""
    ParameterSubset

A struct representing a subset of NMR Hamiltonian parameters to be encoded on qubits.
"""
struct ParameterSubset
    j_couplings::Vector{Tuple{Int,Int}}  # List of (i,j) pairs for J couplings to encode
    h_shifts::Vector{Int}                # List of indices for chemical shifts to encode
end

"""
    QuantumParameterEncoding

A struct representing the encoding of a subset of NMR Hamiltonian parameters onto qubit registers.
"""
struct QuantumParameterEncoding
    j_bits_per_param::Int
    h_bits_per_param::Int
    j_range::Tuple{Float64,Float64}
    h_range::Tuple{Float64,Float64}
    N::Int
    params::ParameterSubset
    
    function QuantumParameterEncoding(
        N::Int,
        params::ParameterSubset;
        j_bits_per_param::Int=4,
        h_bits_per_param::Int=4,
        j_range::Tuple{<:Real,<:Real}=(-1.0, 1.0),
        h_range::Tuple{<:Real,<:Real}=(-1.0, 1.0)
    )
        # Validate parameter indices
        for (i, j) in params.j_couplings
            1 <= i <= N && 1 <= j <= N && i != j || 
                error("Invalid J coupling indices: ($i,$j)")
        end
        
        for i in params.h_shifts
            1 <= i <= N || error("Invalid chemical shift index: $i")
        end
        
        new(j_bits_per_param, h_bits_per_param, 
            Float64.(j_range), Float64.(h_range), N, params)
    end
end

"""
    total_bits(encoding::QuantumParameterEncoding)

Calculate the total number of qubits needed for the parameter encoding.
"""
function total_bits(encoding::QuantumParameterEncoding)
    j_params = length(encoding.params.j_couplings)
    h_params = length(encoding.params.h_shifts)
    
    return (j_params * encoding.j_bits_per_param + 
            h_params * encoding.h_bits_per_param)
end

"""
    binary_to_system(binary_params::Integer, base_system::NMRSystem, encoding::QuantumParameterEncoding)

Apply a binary integer representing parameter subset to create a modified NMRSystem.
"""
function binary_to_system(binary_params::Integer, base_system::NMRSystem, encoding::QuantumParameterEncoding)
    # First, create a copy of the base system
    N = base_system.N
    J = copy(base_system.J)
    h = copy(base_system.h)
    
    j_bits = encoding.j_bits_per_param
    h_bits = encoding.h_bits_per_param
    j_range = encoding.j_range
    h_range = encoding.h_range
    
    # Extract J parameters
    j_mask = (1 << j_bits) - 1
    
    for (idx, (i, j)) in enumerate(encoding.params.j_couplings)
        # Extract j_bits from the binary representation
        j_binary = (binary_params >> ((idx-1) * j_bits)) & j_mask
        
        # Convert to real value
        j_value = binary_to_real(j_binary, j_bits, j_range[1], j_range[2])
        
        # Set the coupling value (ensure symmetric matrix)
        J[i,j] = J[j,i] = j_value
    end
    
    # Extract h parameters
    h_mask = (1 << h_bits) - 1
    h_start_bit = length(encoding.params.j_couplings) * j_bits
    
    for (idx, i) in enumerate(encoding.params.h_shifts)
        # Extract h_bits from the binary representation
        h_binary = (binary_params >> (h_start_bit + (idx-1) * h_bits)) & h_mask
        
        # Convert to real value
        h_value = binary_to_real(h_binary, h_bits, h_range[1], h_range[2])
        
        # Set the chemical shift
        h[i] = h_value
    end
    
    return NMRSystem(N, J, h, base_system.labels, copy(base_system.metadata))
end

"""
    system_to_binary(system::NMRSystem, encoding::QuantumParameterEncoding)

Extract the subset of parameters from an NMRSystem as a binary integer.
"""
function system_to_binary(system::NMRSystem, encoding::QuantumParameterEncoding)
    j_bits = encoding.j_bits_per_param
    h_bits = encoding.h_bits_per_param
    j_range = encoding.j_range
    h_range = encoding.h_range
    
    # Initialize binary parameter
    binary_params = 0
    
    # Encode J parameters
    for (idx, (i, j)) in enumerate(encoding.params.j_couplings)
        # Convert real value to binary
        j_binary = real_to_binary(system.J[i,j], j_bits, j_range[1], j_range[2])
        
        # Add to binary representation
        binary_params |= (j_binary << ((idx-1) * j_bits))
    end
    
    # Encode h parameters
    h_start_bit = length(encoding.params.j_couplings) * j_bits
    
    for (idx, i) in enumerate(encoding.params.h_shifts)
        # Convert real value to binary
        h_binary = real_to_binary(system.h[i], h_bits, h_range[1], h_range[2])
        
        # Add to binary representation
        binary_params |= (h_binary << (h_start_bit + (idx-1) * h_bits))
    end
    
    return binary_params
end

"""
    parameter_index_map(encoding::QuantumParameterEncoding)

Generate a mapping between parameter types/indices and their qubit indices.
"""
function parameter_index_map(encoding::QuantumParameterEncoding)
    j_bits = encoding.j_bits_per_param
    h_bits = encoding.h_bits_per_param
    
    # Initialize map
    param_map = Dict()
    
    # Map J parameters
    bit_offset = 0
    for (idx, (i, j)) in enumerate(encoding.params.j_couplings)
        # For each J parameter, store the qubit indices
        param_map[(:J, i, j)] = collect(bit_offset:(bit_offset + j_bits - 1))
        bit_offset += j_bits
    end
    
    # Map h parameters
    for (idx, i) in enumerate(encoding.params.h_shifts)
        # For each h parameter, store the qubit indices
        param_map[(:h, i, 0)] = collect(bit_offset:(bit_offset + h_bits - 1))
        bit_offset += h_bits
    end
    
    return param_map
end

"""
    generate_all_parameter_configurations(base_system::NMRSystem, encoding::QuantumParameterEncoding)

Generate all possible parameter configurations based on the encoding.
This is for classical simulation/testing of what the quantum superposition represents.
"""
function generate_all_parameter_configurations(base_system::NMRSystem, encoding::QuantumParameterEncoding)
    total_params = 2^total_bits(encoding)
    
    configurations = Vector{NMRSystem}(undef, total_params)
    
    for i in 0:(total_params-1)
        configurations[i+1] = binary_to_system(i, base_system, encoding)
    end
    
    return configurations
end