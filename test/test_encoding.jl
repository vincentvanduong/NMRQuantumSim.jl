using Test
using NMRQuantumSim


@testset "Parameter subset and encoding" begin
    @testset "Encoding" begin
        # Create a base system with 4 spins
        base_system = NMRSystem(4)
        set_coupling!(base_system, 1, 3, 0.3)  # Fixed coupling
        set_coupling!(base_system, 2, 4, -0.2) # Fixed coupling
        set_chemical_shift!(base_system, 1, 0.1) # Fixed shift
        set_chemical_shift!(base_system, 2, 0.2) # Fixed shift
        set_chemical_shift!(base_system, 3, 0.3) # Fixed shift

        # Define which parameters to encode
        params_to_encode = ParameterSubset(
            [(1, 2)],  # Just encode J12
            [4]        # Just encode h4
        )

        # Create the encoding with 3 bits per parameter
        encoding = QuantumParameterEncoding(
            4, # Number of spins
            params_to_encode,
            j_bits_per_param=3, 
            h_bits_per_param=3,
            j_range=(-1.0, 1.0),
            h_range=(-1.0, 1.0)
        )

        # Calculate total qubits needed
        n_qubits = total_bits(encoding)
        @test n_qubits == 6  # Should be 6 (3 for J12 + 3 for h4)

        # Get parameter-to-qubit mapping
        param_map = parameter_index_map(encoding)
        println("J[1,2] is encoded in qubits: ", param_map[(:J, 1, 2)])
        println("h[4] is encoded in qubits: ", param_map[(:h, 4, 0)])

        # Let's test with a specific parameter setting
        test_system = deepcopy(base_system)
        set_coupling!(test_system, 1, 2, 0.5)      # Variable coupling
        set_chemical_shift!(test_system, 4, -0.7)  # Variable shift

        # Convert to binary and back
        binary_params = system_to_binary(test_system, encoding)
        println("Binary representation: ", bitstring(binary_params))

        # Decode and verify
        decoded_system = binary_to_system(binary_params, base_system, encoding)
        println("Original J[1,2]: ", test_system.J[1,2], ", Decoded J[1,2]: ", decoded_system.J[1,2])
        println("Original h[4]: ", test_system.h[4], ", Decoded h[4]: ", decoded_system.h[4])

        @test test_system.J[1,2] ≈ decoded_system.J[1,2] atol=0.15
        @test test_system.h[4] ≈ decoded_system.h[4] atol=0.15

        # Verify fixed parameters were preserved
        println("J[1,3] still fixed at: ", decoded_system.J[1,3])
        println("h[1] still fixed at: ", decoded_system.h[1])

    end
end
