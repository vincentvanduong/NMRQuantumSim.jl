using CairoMakie
using LaTeXStrings
using ColorSchemes

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