using LinearAlgebra
using SparseArrays
using Arpack

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

function heisenberg_hamiltonian(n, J, h)
    # J: n×n coupling matrix where J[i,j] is the coupling between spins i and j
    # h: Array of chemical shifts for each spin
    
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
        H += h[i] * Sz[i]
    end
    
    return H
end

function theta_to_hamiltonian_params(θ, n)
    # θ has length n(n+1)/2 where:
    # - The first n elements are chemical shifts h[1], h[2], ..., h[n]
    # - The remaining n(n-1)/2 elements are the upper triangular elements of J
    
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

function create_heisenberg_hamiltonian(θ, n)
    # Convert θ to J and h
    J, h = theta_to_hamiltonian_params(θ, n)
    
    # Create the Hamiltonian using J and h
    return heisenberg_hamiltonian(n, J, h)
end



function compute_A_mnk(m, n, k, eigenvalues, parameters)
    # Extract parameters
    Δ = parameters.Δ
    l = parameters.l
    
    # Energy difference - transition frequency
    ω_mn = eigenvalues[m] - eigenvalues[n]
    
    # Frequency from phase estimation
    ω_k = 2π * k / (Δ * 2^l)
    
    # Compute A_mnk coefficient using closed-form equation instead of summation
    num = sin(Δ * 2^l * (ω_k - ω_mn) / 2)
    denom = sin(Δ * (ω_k - ω_mn) / 2)

    # Handle the case where denom is very close to zero
    if abs(denom) < 1e-10
        A_mnk = 2^l  # Use the limit value when sin(x)/x approaches 1 as x approaches 0
    else
        prefactor = exp(im * (2.0^l - 1) * Δ * (ω_k - ω_mn) / 2)
        A_mnk = 2.0^(-l) * prefactor * num / denom
    end
    
    return A_mnk
end


function compute_A_tensor_for_parameter(θ, parameters)
    # Compute eigensystem for this parameter set
    eigenvalues, _ = compute_eigensystem(θ, parameters)
    
    # Extract dimensions
    num_states = length(eigenvalues)
    num_k = 2^parameters.l
    
    # Initialize result tensor
    A = Array{ComplexF64}(undef, num_states, num_states, num_k)
    
    # Compute all coefficients
    for m in 1:num_states
        for n in 1:num_states
            for k in 0:(num_k-1)
                A[m, n, k+1] = compute_A_mnk(m, n, k, eigenvalues, parameters)
            end
        end
    end
    
    return A
end

function compute_parallel_A_tensor_for_parameter(θ, parameters)
    return 
end


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

function compute_eigensystem_dense(H::SparseMatrixCSC)
    H_dense = Matrix(H)  # Convert sparse to dense
    eigen_decomp = eigen(Hermitian(H_dense))
    return (eigen_decomp.values, eigen_decomp.vectors)
end


function compute_eigensystem(θ, parameters)
    # Construct Hamiltonian
    H = construct_hamiltonian(θ, parameters)
    
    # Diagonalize using Arpack
    #vals, vecs = eigs(Hermitian(H), nev=min(size(H,1), 10))

    # Diagonalize using Dense
    vals, vecs = compute_eigensystem_dense(H)
    
    return (vals, vecs)
end


function parallel_compute_A_tensors(θ_list, parameters)
    num_θ = length(θ_list)
    
    # First compute one tensor to get dimensions
    test_tensor = compute_A_tensor_for_parameter(θ_list[1], parameters)
    num_states = size(test_tensor, 1)
    num_k = size(test_tensor, 3)
    
    # Initialize 4D result array
    A_tensors = Array{ComplexF64}(undef, num_states, num_states, num_k, num_θ)
    
    # Store the first result
    A_tensors[:, :, :, 1] = test_tensor
    
    return A_tensors
end


# Then handle parallel setup separately
function run_example()
    
    # Define parameters
    parameters = (
        num_spins = 2,      # 1-spin system
        l = 8,              # 6 qubits for phase estimation
        Δ = 2π / 100,       # Time step
    )
    println("Defined parameters")
    
    # Define parameter sets to test
    θ_list = [
        [-0.5, 0.5, 1],  # Example parameters for set 1
        [1.0, 0.5],  # Example parameters for set 2
        [1.0, 1.2]   # Example parameters for set 3
    ]
    
    # Compute A tensors for all parameter sets in parallel
    println("Computing A tensor for first parameter")
    A_tensor = compute_A_tensor_for_parameter(θ_list[1], parameters)
    
    return A_tensor
end

# Execute the example

result = @time run_example()
println(size(result))
println("Computation completed successfully!")