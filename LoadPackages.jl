
using Pkg

# Pkg.add("LinearAlgebra")
# Pkg.add("Random")
# Pkg.add("Statistics")
# Pkg.add("Distributions")
# Pkg.add("InvertedIndices")
# Pkg.add("DataFrames")
# Pkg.add("DataFramesMeta")
# Pkg.add("NamedArrays")
# Pkg.add("CSV")
# Pkg.add("RData")
# Pkg.add("JLD2")
# Pkg.add("ProgressMeter")
# Pkg.add("OnlineStats") # Covariance Matrix Calc
# Pkg.add("MCMCChains")
# Pkg.add("Measures") # for plot specifications
# Pkg.add("StatsPlots")
# Pkg.add("BenchmarkTools")
# Pkg.add("KernelDensity")
# Pkg.add("StatsBase")


using Printf
using LinearAlgebra
using Random
using Statistics
using Distributions
using InvertedIndices
using DataFrames
using DataFramesMeta
using NamedArrays
using CSV
using RData
using JLD2
using ProgressMeter
using OnlineStats # Covariance Matrix Calc
using MCMCChains
using Measures # for plot specifications
using StatsBase
using StatsPlots; gr()
using BenchmarkTools
using KernelDensity
import StatsPlots: cornerplot, cornerplot!


### STRUCTURES

struct Farm_Parish_info
  farm_UID::Int64
  farm_position::Int64
  parish_UID::Int64
  parish_position::Int64
  parish_members_UIDs::Vector{Int64}
  parish_members_positions::Vector{Int64}
end
