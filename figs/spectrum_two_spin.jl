using CairoMakie
using LaTeXStrings
using ColorSchemes

# Define the spectral function
function spectral_function(ω, ω1, ω2, θ, γ)
    # Calculate common term used multiple times
    sqrt_term = sqrt((ω1 - ω2)^2 + 4*θ^2)
    weight1 = 1/4 + θ/(2*sqrt_term)
    weight2 = 1/4 - θ/(2*sqrt_term)
    
    # Term 1
    term1 = (γ * weight1) / (γ^2 + (ω1 + ω2 + 2*θ - sqrt_term - ω)^2)
    
    # Term 2
    term2 = (γ * weight1) / (γ^2 + (ω1 + ω2 - 2*θ + sqrt_term - ω)^2)
    
    # Term 3
    term3 = (γ * weight2) / (γ^2 + (ω1 + ω2 + 2*θ + sqrt_term - ω)^2)
    
    # Term 4
    term4 = (γ * weight2) / (γ^2 + (ω1 + ω2 - 2*θ - sqrt_term + ω)^2)
    
    # Term 5
    term5 = (γ * weight1) / (γ^2 + (ω1 + ω2 + 2*θ - sqrt_term + ω)^2)
    
    # Term 6
    term6 = (γ * weight1) / (γ^2 + (ω1 + ω2 - 2*θ + sqrt_term + ω)^2)
    
    # Term 7
    term7 = (γ * weight2) / (γ^2 + (-ω1 - ω2 + 2*θ + sqrt_term + ω)^2)
    
    # Term 8
    term8 = (γ * weight2) / (γ^2 + (ω1 + ω2 + 2*θ + sqrt_term + ω)^2)
    
    return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8
end

# Set fixed parameters
ω1 = 1.5
ω2 = -0.5
γ = 0.07

# Define a range of θ values to explore
θ_values = [0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2]

# Create frequency range for plotting
ω_range = range(-4.0, 4.0, length=1000)

# Create a figure with appropriate size for publication
fig = Figure(size = (900, 600), fontsize = 16)

# Add a single axis
ax = Axis(
    fig[1, 1],
    xlabel = L"\omega",
    ylabel = L"A(\omega|\theta)",
    title = "Two-spin NMR spectral function",
    titlesize = 20,
    xlabelsize = 18,
    ylabelsize = 18
)

# Select a color palette that works well for scientific publications
colors = ColorSchemes.viridis[range(0.1, 0.9, length=length(θ_values))]

# Plot the spectral function for each θ value
for (i, θ) in enumerate(θ_values)
    spectral_values = [spectral_function(ω, ω1, ω2, θ, γ) for ω in ω_range]
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
    ax, -3.8, 9.0,
    text = L"ω_1 = %$ω1",
    fontsize = 16
)

text!(
    ax, -3.8, 8.5,
    text = L"ω_2 = %$ω2",
    fontsize = 16
)

text!(
    ax, -3.8, 8.0,
    text = L"\gamma = %$γ",
    fontsize = 16
)

# Adjust y-axis limits to accommodate all spectra
ylims!(ax, 0, 10)

# Save the figure in multiple formats
save("figs//two_spin_nmr_spectral_function.eps", fig)
save("figs//two_spin_nmr_spectral_function.pdf", fig)
save("figs//two_spin_nmr_spectral_function.png", fig, px_per_unit = 2)

# Display the figure
fig