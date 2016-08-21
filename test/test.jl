workspace()

include("../src/TraitSimulation.jl")
using DataFrames, TraitSimulation

df = convert(DataFrame, rand(10,3))
names!(df, [:U, :V, :W])
sim_model = Model(T ~ U*V, :CauchitLink, "Normal", [1.0])

simulate(sim_model, df)
