module ExactPottsModel

using Optim
using Combinatorics
using SpecialFunctions
using Random
using JLD2
using Distributed
using Makie
using LaTeXStrings

include("core/NumericalUtils.jl")
include("core/formulas.jl")
include("simulations/Simulations.jl")
include("minimization/minimize_potential.jl")
include("io/jld2_cache.jl")
include("workflows/ConcentrationMinima.jl")
include("workflows/Trajectories.jl")
include("plotting/MakiePlotting.jl")

using .NumericalUtils
using .ExactPottsFormulas
using .Simulations
using .MakiePlotting

end # module
