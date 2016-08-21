workspace()

include("../src/TraitSimulation.jl")
using DataFrames, TraitSimulation

df = convert(DataFrame, rand(10,3))
names!(df, [:U, :V, :W])
sim_model = Model(T ~ U+2.0log(3V+4W)+5.0, "IdentityLink", "Normal", [1.0])

simulate(sim_model, df)
