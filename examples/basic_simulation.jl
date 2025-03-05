using NMRQuantumSim

# Example usage
function run_example()
    # Define parameters for a 3-spin system
    n = 3
    params = create_default_parameters(n)
    
    # Create a coupling matrix (symmetric)
    J = zeros(n, n)
    J[1,2] = J[2,1] = 0.5  # Coupling between spins 1 and 2
    J[1,3] = J[3,1] = 0.3  # Coupling between spins 1 and 3
    J[2,3] = J[3,2] = 0.4  # Coupling between spins 2 and 3
    
    # Chemical shifts for each spin
    h = [2.0, 1.8, 1.6]
    
    # Convert to parameter vector
    θ = vcat(h, [J[i,j] for i in 1:n for j in (i+1):n])
    
    # Compute eigensystem
    eigenvalues, eigenvectors = compute_eigensystem(θ, params)
    println("Eigenvalues: ", eigenvalues)
    
    # Compute A tensor
    A = compute_A_tensor_for_parameter(θ, params)
    println("A tensor shape: ", size(A))
    
    return eigenvalues, A
end

run_example()