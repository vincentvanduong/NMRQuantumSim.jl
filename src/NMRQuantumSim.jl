module NMRQuantumSim

# First include core types that everything else depends on
include("core/types.jl")
include("core/system.jl")

# Then include parameter definitions before they're used
include("simulation/parameters.jl")  # This likely defines ParameterSubset
include("quantum/encoding.jl")       # This defines QuantumParameterEncoding

# Then include calssical inference components
include("classical/posterior_sampling.jl")

# Then include simulation components
include("simulation/operators.jl")
include("simulation/hamiltonian.jl")
include("simulation/evolution.jl")
include("simulation/quantum_fourier_transform.jl")

# Lastly include viualisation components
include("visualisation/posterior_plotting.jl")

# Lastly include IO components
include("io/serialisation.jl")
include("io/data.jl")

# Export public interface
export NMRParameters, create_default_parameters
export ParameterSubset, QuantumParameterEncoding

export BinaryPosteriorSample, convert_binary_sample_to_parameters, collect_posterior_samples
export estimate_posterior_density, marginal_posterior, find_map_estimate, parameter_credible_intervals
export create_spin_operators
export heisenberg_hamiltonian, compute_eigensystem, compute_eigensystem_sparse
export construct_hamiltonian_from_nmr, compute_spectrum, visualise_parameter_space
export compute_qft_array
export compute_A_tensor_for_parameter, parallel_compute_A_tensors

export NMRSystem
export set_coupling!, get_coupling, set_chemical_shift!, get_chemical_shift
export initialize_random_system, hamiltonian_parameters, update_from_parameters!
export binary_to_system, system_to_binary, total_bits, parameter_index_map, real_to_binary, binary_to_real
export save_system, load_system, calculate_spectrum

export plot_spectral_map
export save_experiment, load_experiment, list_experiments, run_and_save_experiment
export plot_1d_density
export create_corner_plot

end # module