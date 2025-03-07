using CairoMakie
using LaTeXStrings
using ColorSchemes

function plot_spectral_map(spectrum_2d, parameters, file_name)

    # Extract parameters
    #Δ = parameters.Δ
    #l = parameters.l

    n_freq, n_configs = size(spectrum_2d)

    freq_range = range(0, n_freq-1, length=n_freq)
    θ_range = range(0, n_configs-1, length=n_configs)

    # Create the figure
    fig = Figure(size = (900, 600), fontsize = 16)

    # Create a single axis for the heatmap
    ax = Axis(
        fig[1, 1],
        xlabel = L"k",
        ylabel = L"\theta",
        title = "Parameter-dependent NMR spectral map",
        titlesize = 20,
        xlabelsize = 18,
        ylabelsize = 18
    )

    # Create heatmap of the spectral map
    hm = heatmap!(
        ax, 
        θ_range,
        freq_range,
        transpose(spectrum_2d)  # Transpose to match the desired orientation
    )

    hm.colormap = :viridis

    # Add colorbar
    colorbar = Colorbar(
        fig[1, 2],
        hm,
        label = L"\mathrm{Spectral \; Intensity}",
        labelsize = 16,
        ticklabelsize = 14
    )

    # Add parameter annotations
    #text!(
    #    ax, -3.8, 1.4,
    #    text = L"h_1 = %$h1",
    #    fontsize = 16
    #)

    #text!(
    #    ax, -3.8, 1.3,
    #    text = L"h_2 = %$h2",
    #    fontsize = 16
    #)

    #text!(
    #    ax, -3.8, 1.2,
    #    text = L"\gamma = %$γ",
    #    fontsize = 16
    #)

    # Save the figure
    base_dir = "figs"

    # Create the full paths for different formats
    eps_path = joinpath(base_dir, file_name * ".eps")
    pdf_path = joinpath(base_dir, file_name * ".pdf")
    png_path = joinpath(base_dir, file_name * ".png")

    save(eps_path, fig)
    save(pdf_path, fig)
    save(png_path, fig, px_per_unit = 2)

    # Display the figure
    fig
end


"""
    plot_parameter_posterior(posterior_samples; title="Parameter Posterior", n_bins=30, kwargs...)

Plot the posterior distribution of parameters.

# Arguments
- `posterior_samples::Dict`: Dictionary with parameter names as keys and sample vectors as values
- `title::String`: Plot title (default: "Parameter Posterior")
- `n_bins::Int`: Number of bins for histogram (default: 30)
- `kwargs...`: Additional keyword arguments for the plot

# Returns
- A Plots.jl plot object
"""
function plot_parameter_posterior(posterior_samples; title="Parameter Posterior", n_bins=30, kwargs...)
    n_params = length(posterior_samples)
    
    if n_params > 9
        println("Warning: Plotting many parameters. Consider selecting fewer parameters.")
    end
    
    # Calculate layout rows and columns
    n_cols = min(3, n_params)
    n_rows = ceil(Int, n_params / n_cols)
    
    p = plot(layout=(n_rows, n_cols), size=(300*n_cols, 250*n_rows), title=title)
    
    for (i, (param_name, samples)) in enumerate(posterior_samples)
        histogram!(p[i], samples, 
            bins=n_bins,
            normalize=true,
            xlabel=param_name,
            ylabel= i % n_cols == 1 ? "Density" : "",
            title="",
            label=false,
            fillalpha=0.7,
            titlefontsize=10,
            kwargs...
        )
        
        # Add mean and 95% CI
        mean_val = mean(samples)
        ci_low, ci_high = quantile(samples, [0.025, 0.975])
        
        vline!(p[i], [mean_val], linestyle=:dash, linewidth=2, color=:black, label="")
        vspan!(p[i], [ci_low, ci_high], alpha=0.2, color=:blue, label="")
    end
    
    return p
end