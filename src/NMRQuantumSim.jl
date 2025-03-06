module NMRQuantumSim

# Include all component files
include("core/types.jl")
include("core/system.jl")
# include("core/utils.jl")
include("simulation/parameters.jl")
include("simulation/operators.jl")
include("simulation/hamiltonian.jl")
include("simulation/evolution.jl")
include("simulation/quantum_fourier_transform.jl")
include("quantum/encoding.jl")
include("quantum/circuits.jl")
include("quantum/measurement.jl")
include("classical/spectrum.jl")
#include("classical/inference.jl")
include("io/serialisation.jl")
include("io/visualisation.jl")

# Export public interface
export NMRParameters, create_default_parameters
export create_spin_operators
export heisenberg_hamiltonian, compute_eigensystem, compute_eigensystem_sparse
export construct_hamiltonian_from_nmr, compute_qft_array
export compute_A_tensor_for_parameter, parallel_compute_A_tensors
export NMRSystem
export set_coupling!, get_coupling, set_chemical_shift!, get_chemical_shift
export initialize_random_system, hamiltonian_parameters, update_from_parameters!
export QuantumParameterEncoding, ParameterSubset
export binary_to_system, system_to_binary, total_bits, parameter_index_map, real_to_binary, binary_to_real
export save_system, load_system, calculate_spectrum
export plot_spectral_map

end # module