workspace()

include("../src/TraitSimulation.jl")

using DataFrames, TraitSimulation

# create data frame for testing
df = convert(DataFrame, rand(10,6))
names!(df, [:A, :B, :C, :D, :E, :F])

# simulate from the model
formula = T ~ A+2B*C+3log(D*(E+F))
sim_model = Model(formula, IdentityLink(), TResponse(1))
simulate(sim_model, df)
