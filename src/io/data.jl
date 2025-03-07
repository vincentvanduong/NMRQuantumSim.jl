# File: nmr_data_framework.jl

using HDF5
using Plots
using Statistics
using Random
using LinearAlgebra
using DataFrames
using StatsBase
import Dates

"""
    save_experiment(name, parameters, results; path="./data", metadata=Dict())

Save experimental results to an HDF5 file with proper naming and organization.

# Arguments
- `name::String`: Name of the experiment
- `parameters::Dict`: Dictionary of parameters used in the experiment
- `results::Dict`: Dictionary of result data to save
- `path::String`: Directory path to save the file (default: "./data")
- `metadata::Dict`: Additional metadata about the experiment (optional)

# Returns
- `filename::String`: Full path of the saved file
"""
function save_experiment(name, parameters, results; path="./data", metadata=Dict())
    # Create directory if it doesn't exist
    mkpath(path)
    
    # Create a timestamp
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
    
    # Create a unique filename
    filename = joinpath(path, "$(name)_$(timestamp).h5")
    
    # Save to HDF5
    h5open(filename, "w") do file
        # Create a group for parameters
        params_group = create_group(file, "parameters")
        for (key, value) in parameters
            params_group[key] = value
        end
        
        # Create a group for results
        results_group = create_group(file, "results")
        for (key, value) in results
            results_group[key] = value
        end
        
        # Add metadata
        meta_group = create_group(file, "metadata")
        
        # Add standard metadata
        meta_group["timestamp"] = timestamp
        meta_group["experiment_name"] = name
        
        # Add custom metadata
        for (key, value) in metadata
            meta_group[key] = value
        end
    end
    
    println("Experiment saved to $filename")
    return filename
end

"""
    load_experiment(filename)

Load an experiment from an HDF5 file.

# Arguments
- `filename::String`: Path to the HDF5 file

# Returns
- `experiment::Dict`: Dictionary containing parameters, results, and metadata
"""
function load_experiment(filename)
    experiment = Dict(
        "parameters" => Dict(),
        "results" => Dict(),
        "metadata" => Dict()
    )
    
    h5open(filename, "r") do file
        # Load parameters
        params_group = file["parameters"]
        for key in keys(params_group)
            experiment["parameters"][key] = read(params_group[key])
        end
        
        # Load results
        results_group = file["results"]
        for key in keys(results_group)
            experiment["results"][key] = read(results_group[key])
        end
        
        # Load metadata
        meta_group = file["metadata"]
        for key in keys(meta_group)
            experiment["metadata"][key] = read(meta_group[key])
        end
    end
    
    return experiment
end

"""
    list_experiments(path="./data")

List all experiments in the given directory.

# Arguments
- `path::String`: Directory path to search (default: "./data")

# Returns
- `experiments::DataFrame`: DataFrame with experiment details
"""
function list_experiments(path="./data")
    files = filter(f -> endswith(f, ".h5"), readdir(path, join=true))
    
    if isempty(files)
        println("No experiments found in $path")
        return DataFrame()
    end
    
    # Prepare data for DataFrame
    names = String[]
    timestamps = String[]
    parameter_counts = Int[]
    result_keys = Vector{Vector{String}}()
    
    for file in files
        experiment = load_experiment(file)
        push!(names, experiment["metadata"]["experiment_name"])
        push!(timestamps, experiment["metadata"]["timestamp"])
        push!(parameter_counts, length(keys(experiment["parameters"])))
        push!(result_keys, collect(keys(experiment["results"])))
    end
    
    return DataFrame(
        Filename = files,
        Name = names,
        Timestamp = timestamps,
        Parameters = parameter_counts,
        Results = result_keys
    )
end

"""
    run_and_save_experiment(name, experiment_fn, parameters; path="./data", metadata=Dict())

Run an experiment function with the given parameters and save the results.

# Arguments
- `name::String`: Name of the experiment
- `experiment_fn::Function`: Function that takes parameters and returns results dict
- `parameters::Dict`: Dictionary of parameters to pass to the experiment function
- `path::String`: Directory path to save the file (default: "./data")
- `metadata::Dict`: Additional metadata about the experiment (optional)

# Returns
- `filename::String`: Full path of the saved file
- `results::Dict`: Dictionary of results from the experiment
"""
function run_and_save_experiment(name, experiment_fn, parameters; path="./data", metadata=Dict())
    # Run the experiment
    println("Running experiment: $name")
    results = experiment_fn(parameters)
    
    # Save the experiment
    filename = save_experiment(name, parameters, results, path=path, metadata=metadata)
    
    return filename, results
end