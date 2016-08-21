workspace()

include("../src/TraitSimulation.jl")

using DataFrames, TraitSimulation

# create data frame for testing
df = convert(DataFrame, rand(10,3))
names!(df, [:U, :V, :W])

# simulate from the model
normal_dist = Dict(:Name => :Normal, :σ => 1.0)
sim_model = Model(T ~ U*2V+2log(W), :CauchitLink, normal_dist)
simulate(sim_model, df)
