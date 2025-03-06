# test_types.jl
module TestTypes

struct ParameterSubset
    j_couplings::Vector{Tuple{Int,Int}}
    h_shifts::Vector{Int}
end

struct QuantumParameterEncoding
    j_bits_per_param::Int
    h_bits_per_param::Int
    j_range::Tuple{Float64,Float64}
    h_range::Tuple{Float64,Float64}
    N::Int
    params::ParameterSubset
    
    function QuantumParameterEncoding(
        N::Int,
        params::ParameterSubset;
        j_bits_per_param::Int=4,
        h_bits_per_param::Int=4,
        j_range::Tuple{<:Real,<:Real}=(-1.0, 1.0),
        h_range::Tuple{<:Real,<:Real}=(-1.0, 1.0)
    )
        new(j_bits_per_param, h_bits_per_param, 
            Float64.(j_range), Float64.(h_range), N, params)
    end
end

function test_create()
    params = ParameterSubset([(1,2)], [3])
    encoding = QuantumParameterEncoding(3, params)
    println("Successfully created encoding")
    return encoding
end

end # module

# Test outside the module
using .TestTypes
TestTypes.test_create()