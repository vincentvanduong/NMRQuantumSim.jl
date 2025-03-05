using CairoMakie
using LaTeXStrings
using ColorSchemes

# Define the spectral function
function spectral_function(ω, ω0, θ, γ)
    term1 = 2 * (ω0^2)/(θ^2 + ω0^2) * (γ/(ω^2 + γ^2))
    term2 = (θ^2)/(θ^2 + ω0^2) * (
        γ/((ω + 2 * sqrt(ω0^2 + θ^2))^2 + γ^2) + 
        γ/((ω - 2 * sqrt(ω0^2 + θ^2))^2 + γ^2)
    )
    return (term1 + term2)/2
end

# Set fixed parameters
ω0 = 1.0    # central frequency
γ = 0.05    # damping/linewidth

# Define a range of θ values to explore
θ_values = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5]

# Create frequency range for plotting
ω_range = range(-4.0, 4.0, length=1000)

# Create a figure with appropriate size for publication
fig = Figure(size = (800, 600), fontsize = 16)

# Add a single axis
ax = Axis(
    fig[1, 1],
    xlabel = L"\omega",
    ylabel = L"A(\omega|\theta)",
    title = "NMR spectral functions for different coupling strengths",
    titlesize = 20,
    xlabelsize = 18,
    ylabelsize = 18
)

# Select a color palette that works well for scientific publications
colors = ColorSchemes.viridis[range(0.1, 0.9, length=length(θ_values))]

# Plot the spectral function for each θ value
for (i, θ) in enumerate(θ_values)
    spectral_values = [spectral_function(ω, ω0, θ, γ) for ω in ω_range]
    lines!(ax, ω_range, spectral_values, 
          color = colors[i], 
          linewidth = 2,
          label = L"\theta = %$(round(θ, digits=1))")
end

# Add a legend with correct position syntax
axislegend(ax, position = :rt, nbanks = 2)  # rt means right-top, nbanks creates two columns

# Add grid lines for better readability
ax.xgridvisible = false
ax.ygridvisible = false

# Add annotations about fixed parameters
text!(
    ax, -3.9, 0.9*25,
    text = L"\omega_0 = %$ω0",
    fontsize = 16
)

text!(
    ax, -3.9, 0.85*25,
    text = L"\gamma = %$γ",
    fontsize = 16
)

# Adjust y-axis limits to accommodate all spectra
ylims!(ax, 0, 25)

# Save the figure in multiple formats
save("figs//nmr_spectral_function_family.pdf", fig)
save("figs//nmr_spectral_function_family.eps", fig)
save("figs//nmr_spectral_function_family.png", fig, px_per_unit = 2)

# Display the figure
fig