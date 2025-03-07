using NMRQuantumSim

# List all experiments
experiments_df = list_experiments()
println(experiments_df)

# Load a specific experiment and analyze
loaded_experiment = load_experiment("data/two_spin_simulation_2025-03-06_21-48-08.h5")
println(loaded_experiment)

parameters = loaded_experiment["results"]["parameters"]
posterior = loaded_experiment["results"]["posterior"]

fig = plot_1d_density(parameters, posterior)

## Save the figure in multiple formats
#save("figs/two_spin_simulation_posterior.eps", fig)
#save("figs/two_spin_simulation_posterior.pdf", fig)
#save("figs/two_spin_simulation_posterior.png", fig, px_per_unit = 2)

#fig