workspace()

include("../src/TraitSimulation.jl")

using DataFrames, TraitSimulation

# create data frame for testing
df = convert(DataFrame, rand(10,3))
names!(df, [:U, :V, :W])

# simulate from the model
sim_model = Model(T ~ U*2V+2log(W), :CauchitLink, NormalResponse(1.0))
simulate(sim_model, df)
