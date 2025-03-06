using LinearAlgebra


function construct_hamiltonian_from_nmr(system::NMRSystem)
    N = system.N
    dim = 2^N

    # Create spin operators
    Sx, Sy, Sz = create_spin_operators(N)
    
    # Initialize Hamiltonian matrix
    H = spzeros(ComplexF64, dim, dim)
    
    # For each spin pair (i,j), add J_ij * S_i · S_j term
    for i in 1:N
        for j in (i+1):N
            J_ij = system.J[i,j]
            
            # Skip if coupling is zero
            isapprox(J_ij, 0.0, atol=1e-10) && continue
            
            H += J_ij * (Sx[i]*Sx[j] + Sy[i]*Sy[j] + Sz[i]*Sz[j])
        end
    end
    
    # Add chemical shift terms: ∑_i h_i S_i^x
    for i in 1:N
        h_i = system.h[i]
        
        # Skip if chemical shift is zero
        isapprox(h_i, 0.0, atol=1e-10) && continue
                
        # Add to Hamiltonian
        H += h_i * Sx[i]
    end
    
    return H
end