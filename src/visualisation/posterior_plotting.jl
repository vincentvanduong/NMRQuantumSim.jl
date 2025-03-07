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

# Define the spectral function
function plot_1d_density(x, y)
    
    # Create a figure with appropriate size for publication
    fig = Figure(size = (900, 600), fontsize = 16)

    # Add a single axis
    ax = Axis(
        fig[1, 1],
        xlabel = L"\theta",
        ylabel = L"P(\theta|X)",
        title = "Two-spin NMR posterior",
        titlesize = 20,
        xlabelsize = 18,
        ylabelsize = 18
    )

    # Select a color palette that works well for scientific publications
    #colors = ColorSchemes.viridis[range(0.1, 0.9, length=length(θ_values))]

    barplot!(ax, y, x,
            gap = 0,
            color = :gray85,
            strokecolor = :black,
            strokewidth = 1)


    #    line!color = colors[i], 
        linewidth = 2
        #label = L"\theta = %$(round(θ, digits=1))")
    #)

    # Add a legend with correct position syntax
    #axislegend(ax, position = :rt, nbanks = 2)  # rt means right-top, nbanks creates two columns

    # Add grid lines for better readability
    ax.xgridvisible = false
    ax.ygridvisible = false

    # Add annotations about fixed parameters
    #text!(
    #    ax, -3.8, 9.0,
    #    text = L"ω_1 = %$ω1",
    #    fontsize = 16
    #)

    #text!(
    #    ax, -3.8, 8.5,
    #    text = L"ω_2 = %$ω2",
    #    fontsize = 16
    #)

    #text!(
    #    ax, -3.8, 8.0,
    #    text = L"\gamma = %$γ",
    #    fontsize = 16
    #)

    # Adjust y-axis limits to accommodate all spectra
    #ylims!(ax, 0, 10)

    return fig
end