#module TraitSimulationTest

include("../src/TraitSimulation.jl")

using DataFrames, TraitSimulation

if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
end

# create data frame for testing
srand(1)
data = rand(5,6)
df = convert(DataFrame, data)
cov_mat = cor(data')
names!(df, [:A, :B, :C, :D, :E, :F])

# simulate a single trait using GLM
sim_model = Model(T ~ A + 2B*C, IdentityLink(), NormalResponse(1.0))
y = simulate(sim_model, df)

# simulate two traits with the same link and response using GLM
formulae = [T1 ~ A+2B*C, T2 ~ A + 2log(B+C)+2.0]
sim_model = Model(formulae, IdentityLink(), NormalResponse(1.0))
y = simulate(sim_model, df)

# simulate three traits with different link and response using GLM
formulae = [T1 ~ A+2B*C, T2 ~ A+2log(B+C)+2.0, T3 ~ A+B+C+1.0]
links = [IdentityLink(), LogitLink(), LogLink()]
resp_dists = [NormalResponse(1.0), BinomialResponse(100), PoissonResponse()]
sim_model = Model(formulae, links, resp_dists)
y = simulate(sim_model, df)

# simulate a single trait with two variance components using GLMM
vc = [VarianceComponent(0.2, cov_mat), VarianceComponent(0.8, eye(5))]
sim_model = Model(T ~ A+2B*C, LogitLink(), BinomialResponse(100), vc)
y = simulate(sim_model, df)

# simulate two traits with two variance components using GLMM
# the two traits can have different response distribution
vc = [VarianceComponent([0.2, 0.3], cov_mat),
      VarianceComponent([0.8, 0.7], eye(5))]
formulae = [T1 ~ A+2B*C, T2 ~ C+log(C)+3.0]
links = [IdentityLink(), LogitLink()]
resp_dists = [NormalResponse(1.0), PoissonResponse()]
sim_model = Model(formulae, links, resp_dists, vc)
y = simulate(sim_model, df)

# simulate two traits with two variance components with cross covariances
# using GLMM
vc = [VarianceComponent([0.2 -0.1; -0.1 0.3], cov_mat),
      VarianceComponent([0.8 -0.2; -0.2 0.7], eye(5))]
formulae = [T1 ~ A+2B*C, T2 ~ C+log(C)+3.0]
sim_model = Model(formulae, IdentityLink(), NormalResponse(1.0), vc)
y = simulate(sim_model, df)
println(y)

#end
