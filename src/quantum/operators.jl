using SparseArrays, LinearAlgebra

"""
    create_spin_operators(n)

Create spin-1/2 operators for an n-spin system.
Returns dictionaries of Sx, Sy, Sz operators for each site.
"""
function create_spin_operators(n)
    # Single-site operators
    σx = sparse([0 1; 1 0])
    σy = sparse([0 -im; im 0])
    σz = sparse([1 0; 0 -1])
    I2 = sparse(I, 2, 2)
    
    # Initialize empty dictionaries for operators
    Sx, Sy, Sz = Dict(), Dict(), Dict()
    
    # Construct operators for each site
    for i in 1:n
        # Create identity operators for all sites except i
        terms = [j == i ? σx : I2 for j in 1:n]
        Sx[i] = kron(terms...)
        
        terms = [j == i ? σy : I2 for j in 1:n]
        Sy[i] = kron(terms...)
        
        terms = [j == i ? σz : I2 for j in 1:n]
        Sz[i] = kron(terms...)
    end
    
    return Sx, Sy, Sz
end