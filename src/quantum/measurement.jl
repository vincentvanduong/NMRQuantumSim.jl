"""
    decode_binary_measurement(binary_result::Integer, base_system::NMRSystem, encoding::QuantumParameterEncoding)

Decode a measurement result from the quantum register into actual Hamiltonian parameters.
"""
function decode_binary_measurement(binary_result::Integer, base_system::NMRSystem, encoding::QuantumParameterEncoding)
    return binary_to_system(binary_result, base_system, encoding)
end

# Other measurement-related functions...