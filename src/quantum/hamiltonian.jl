using SparseArrays, LinearAlgebra, Arpack

"""
    heisenberg_hamiltonian(n, J, h)

Construct the Heisenberg Hamiltonian for an n-spin system
with coupling matrix J and chemical shifts h.
"""
function heisenberg_hamiltonian(n, J, h)
    # Create spin operators
    Sx, Sy, Sz = create_spin_operators(n)
    
    # Initialize Hamiltonian
    dim = 2^n
    H = spzeros(ComplexF64, dim, dim)
    
    # Add coupling terms
    for i in 1:n
        for j in i+1:n
            J_ij = J[i,j]
            H += J_ij * (Sx[i]*Sx[j] + Sy[i]*Sy[j] + Sz[i]*Sz[j])
        end
    end
    
    # Add chemical shift terms
    for i in 1:n
        H += h[i] * Sx[i]
    end
    
    return H
end

"""
    theta_to_hamiltonian_params(θ, n)

Convert a parameter vector θ to coupling matrix J and chemical shifts h.
"""
function theta_to_hamiltonian_params(θ, n)
    # Extract chemical shifts
    h = θ[1:n]
    
    # Create coupling matrix
    J = zeros(n, n)
    
    # Fill coupling matrix using the remaining parameters
    idx = n + 1
    for i in 1:n
        for j in i+1:n
            J[i,j] = J[j,i] = θ[idx]
            idx += 1
        end
    end
    
    return J, h
end

"""
    create_heisenberg_hamiltonian(θ, n)

Create a Heisenberg Hamiltonian from parameter vector θ for n spins.
"""
function create_heisenberg_hamiltonian(θ, n)
    # Convert θ to J and h
    J, h = theta_to_hamiltonian_params(θ, n)
    
    # Create the Hamiltonian using J and h
    return heisenberg_hamiltonian(n, J, h)
end

"""
    compute_eigensystem(H::SparseMatrixCSC)

Compute the eigenvalues and eigenvectors of the sparse Hamiltonian.
"""
function compute_eigensystem(H::SparseMatrixCSC)
    # For large matrices, just compute a subset of eigenvalues/vectors
    # Change nev parameter based on how many eigenvalues you need
    vals, vecs = eigs(Hermitian(H), nev=min(size(H,1), 100))
    return (vals, vecs)
end

"""
    compute_eigensystem(H::AbstractMatrix)

Compute the eigenvalues and eigenvectors of a dense Hamiltonian.
"""
function compute_eigensystem(H::AbstractMatrix)
    eigen_decomp = eigen(Hermitian(H))
    vals, vecs = eigen_decomp.values, eigen_decomp.vectors
    return (vals, vecs)
end

"""
    construct_hamiltonian_from_system(system::NMRSystem)

Construct the Hamiltonian matrix from an NMR system with couplings and chemical shifts.
"""
function construct_hamiltonian_from_system(system::NMRSystem)
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