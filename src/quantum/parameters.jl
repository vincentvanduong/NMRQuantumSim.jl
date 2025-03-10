"""
    NMRParameters

A struct containing simulation parameters for NMR quantum computation.
"""
struct NMRParameters
    num_spins::Int               # Number of nuclear spins
    l::Int                       # Number of phase estimation qubits
    Δ::Float64                   # Time evolution step
    num_m::Int                   # Dimension of first register
    num_n::Int                   # Dimension of second register
end

"""
    create_default_parameters(num_spins)

Create a default parameter set for an NMR system with num_spins nuclear spins.
"""
function create_default_parameters(num_spins)
    l = 6  # 10 qubits for phase estimation
    Δ = 2π / 1000  # Time step
    dim = 2^num_spins  # Hilbert space dimension
    
    return NMRParameters(
        num_spins,
        l,
        Δ,
        dim,
        dim
    )
end