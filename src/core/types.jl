using StaticArrays
using LinearAlgebra

"""
    AbstractNMRSystem

Abstract type for NMR spin systems.
"""
abstract type AbstractNMRSystem end

"""
    ParameterSubset

A struct representing a subset of NMR Hamiltonian parameters.
"""
struct ParameterSubset
    j_couplings::Vector{Tuple{Int,Int}}  # List of (i,j) pairs for J couplings
    h_shifts::Vector{Int}                # List of indices for chemical shifts
end