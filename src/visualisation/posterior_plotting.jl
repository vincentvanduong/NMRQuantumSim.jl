using CairoMakie, Statistics

"""
    create_corner_plot(samples::Vector{BinaryPosteriorSample}, param_names::Vector{Symbol}; 
                      fig_size=(900, 900), point_size=5, color_scale=:viridis)

Create a corner plot (pair plot) showing the 2D marginal distributions of parameters.
Uses raw log probabilities for coloring without smoothing.

Returns a Figure object that can be saved or displayed.
"""
function create_corner_plot(samples::Vector{BinaryPosteriorSample}, param_names::Vector{Symbol}; 
                          fig_size=(900, 900), point_size=5, color_scale=:viridis)
    # Extract parameter values and log probabilities
    n_params = length(param_names)
    n_samples = length(samples)
    
    # Create a matrix of parameter values
    param_values = zeros(n_samples, n_params)
    for i in 1:n_samples
        for j in 1:n_params
            param_values[i, j] = samples[i].parameter_values[param_names[j]]
        end
    end
    
    # Get log probabilities for coloring
    log_probs = [s.log_probability for s in samples]
    
    # Create figure with grid of plots
    fig = CairoMakie.Figure(size=fig_size)
    
    # Find parameter ranges for consistent axes
    param_ranges = [(minimum(param_values[:, i]), maximum(param_values[:, i])) for i in 1:n_params]
    
    # Create corner plot layout
    for i in 1:n_params
        for j in 1:n_params
            # Only create plots for lower triangular part (including diagonal)
            if j > i
                continue
            end
            
            # Create axis
            ax = CairoMakie.Axis(fig[i, j])
            
            if i == j
                # Diagonal: 1D histogram
                CairoMakie.hist!(ax, param_values[:, i], bins=20)
                ax.title = string(param_names[i])
            else
                # Off-diagonal: 2D scatter plot
                CairoMakie.scatter!(ax, param_values[:, j], param_values[:, i], 
                         color=log_probs, colormap=color_scale, 
                         markersize=point_size)
                
                # Set axis limits
                CairoMakie.xlims!(ax, param_ranges[j])
                CairoMakie.ylims!(ax, param_ranges[i])
            end
            
            # Only show labels on the edge plots
            if i < n_params
                CairoMakie.hidexdecorations!(ax, grid=false, ticks=false)
            else
                ax.xlabel = string(param_names[j])
            end
            
            if j > 1
                CairoMakie.hideydecorations!(ax, grid=false, ticks=false)
            else
                ax.ylabel = string(param_names[i])
            end
        end
    end
    
    # Add colorbar for log probability
    CairoMakie.Colorbar(fig[1:n_params, n_params+1], colormap=color_scale, 
             label="Log Probability", limits=(minimum(log_probs), maximum(log_probs)))
    
    return fig
end