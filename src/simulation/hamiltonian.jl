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
function theta_to_hamiltonian_params(θv, n)
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
    construct_hamiltonian(θ, parameters)

Construct the Hamiltonian for simulation parameters.
"""
function construct_hamiltonian(θ, parameters)
    # Extract system size
    num_spins = parameters.num_spins
    
    H = create_heisenberg_hamiltonian(θ, num_spins)
    
    return H
end


function compute_eigensystem_sparse(H::SparseMatrixCSC)
    # For large matrices, just compute a subset of eigenvalues/vectors
    # Change nev parameter based on how many eigenvalues you need
    vals, vecs = eigs(Hermitian(H), nev=min(size(H,1), 10))
    return (vals, vecs)
end

"""
    compute_eigensystem_dens(H::SparseMatrixCSC)

Compute the eigenvalues and eigenvectors of the sparse Hamiltonian.
"""
function compute_eigensystem_dense(H::SparseMatrixCSC)
    H_dense = Matrix(H)  # Convert sparse to dense
    eigen_decomp = eigen(Hermitian(H_dense))
    vals, vecs = eigen_decomp.values, eigen_decomp.vectors
    return (vals, vecs)
end

"""
    compute_eigensystem(θ, parameters)

Compute the eigenvalues and eigenvectors of the Hamiltonian.
"""
function compute_eigensystem(θ, parameters)
    # Construct sparse Hamiltonian
    H = construct_hamiltonian(θ, parameters)
    
    # Diagonalise a dense Hamiltonian (not using Arpack)
    vals, vecs = compute_eigensystem_dense(H)
    
    return (vals, vecs)
end
