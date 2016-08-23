workspace()

include("../src/TraitSimulation.jl")

using DataFrames, TraitSimulation

# create data frame for testing
df = convert(DataFrame, rand(10,6))
names!(df, [:A, :B, :C, :D, :E, :F])

# simulate from the model
formulae = [T1 ~ A+2B*C+3log(D*(E+F)), T2 ~ C + 2D]
links = [IdentityLink(), LogitLink()]
resp_dist = [NormalResponse(1.0), BinomialResponse(100)]
sim_model = Model(formulae, links, resp_dist)
y = simulate(sim_model, df)
println(y)
