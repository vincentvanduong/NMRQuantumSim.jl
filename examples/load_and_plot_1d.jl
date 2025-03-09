using NMRQuantumSim
using CairoMakie

# List all experiments
experiments_df = list_experiments()
println(experiments_df)

# Load a specific experiment and analyze
loaded_experiment = load_experiment("data/two_spin_simulation_2025-03-07_23-13-08.h5")
#println(loaded_experiment)

binary_samples = loaded_experiment["results"]["binary_samples"]
posterior_prob = loaded_experiment["results"]["posterior_prob"]

fig = plot_1d_density(binary_samples, posterior_prob)

## Save the figure in multiple formats
#CairoMakie.save("figs/two_spin_simulation_posterior.eps", fig)
#CairoMakie.save("figs/two_spin_simulation_posterior.pdf", fig)
#CairoMakie.save("figs/two_spin_simulation_posterior.png", fig, px_per_unit = 2)

#fig

#binary_samples