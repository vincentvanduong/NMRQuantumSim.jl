module NMRQuantumSim

# First include core types that everything else depends on
include("core/types.jl")
include("core/system.jl")

# Then include parameter definitions before they're used
include("quantum/parameters.jl")
include("quantum/encoding.jl")
include("quantum/operators.jl")
include("quantum/hamiltonian.jl")
include("quantum/measurement.jl")
include("quantum/circuits.jl")

# Then include calssical inference components
include("inference/posterior_sampling.jl")

# Then include simulation components
include("simulation/evolution.jl")
include("simulation/phase_estimation.jl")

# Then include io components
include("io/serialisation.jl")

# Export public interface
export NMRSystem

export NMRParameters, create_default_parameters
export ParameterSubset, QuantumParameterEncoding
export binary_to_system, system_to_binary, total_bits, parameter_index_map, real_to_binary, binary_to_real

export set_coupling!, get_coupling, set_chemical_shift!, get_chemical_shift
export initialize_random_system, hamiltonian_parameters, update_from_parameters!

export heisenberg_hamiltonian, compute_eigensystem, construct_hamiltonian_from_system, create_spin_operators

export compute_qft_array, compute_spectrum, visualise_parameter_space

export BinaryPosteriorSample, convert_binary_sample_to_parameters, collect_posterior_samples

end # module