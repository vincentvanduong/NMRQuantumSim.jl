using NMRQuantumSim
using CairoMakie

# List all experiments
experiments_df = list_experiments()
println(experiments_df)

# Load a specific experiment and analyze
loaded_experiment = load_experiment("data/two_spin_simulation_2025-03-07_23-16-47.h5")

# Create a base system with 2 spins
N_spins = loaded_experiment["parameters"]["N_spins"]
chemical_shift1 = loaded_experiment["parameters"]["chemical_shift1"]
chemical_shift2 = loaded_experiment["parameters"]["chemical_shift2"]

# Create a base system with 4 spins
base_system = NMRSystem(N_spins)
set_chemical_shift!(base_system, 1, chemical_shift1) # Fixed shift
set_chemical_shift!(base_system, 2, chemical_shift2) # Fixed shift

# Define parameters to encode (J12 and h4)
params_to_encode = ParameterSubset(
    [(1, 2)],  # Just encode J12
    [1]        # Just encode h1
)

# Create the encoding with 3 bits per parameter
encoding = QuantumParameterEncoding(
    N_spins, 
    params_to_encode,
    j_bits_per_param=8, 
    h_bits_per_param=8,
    j_range=(0, 1),
    h_range=(-1, 1)
)

# Assume experiement contains binary samples and their log probabilities
binary_samples = Int64.([param for param in loaded_experiment["results"]["binary_samples"]])
log_probs = [log(prob) for prob in loaded_experiment["results"]["posterior_prob"]]

# Convert samples to parameter space
posterior_samples = collect_posterior_samples(
    binary_samples, log_probs, base_system, encoding
)

# Create list of parameter names we're interested in
param_names = [Symbol("J_1_2"), Symbol("h_1")]

# Create and display the corner plot
fig = create_corner_plot(posterior_samples, param_names)
fig