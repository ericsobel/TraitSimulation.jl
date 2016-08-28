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
GRM = cor(data')
names!(df, [:A, :B, :C, :D, :E, :F])

# simulate a single trait using GLM
sim_model = Model(T ~ A+2B*C, IdentityLink(), NormalResponse(1.0))
y = simulate(sim_model, df)

# simulate two traits with the same link and response using GLM
μ = [T1 ~ A+2B*C, T2 ~ A+2log(B+C)+2.0]
sim_model = Model(μ, IdentityLink(), NormalResponse(1.0))
y = simulate(sim_model, df)

# simulate three traits with different link and response using GLM
μ = [T1 ~ A+2B*C, T2 ~ A+2log(B+C)+2.0, T3 ~ A+B+C+1.0]
links = [IdentityLink(), LogitLink(), LogLink()]
resp_dists = [NormalResponse(1.0), BinomialResponse(100), PoissonResponse()]
sim_model = Model(μ, links, resp_dists)
y = simulate(sim_model, df)

# simulate a single trait with two variance components using GLMM
Σ = [VarianceComponent(0.2, GRM), VarianceComponent(0.8, eye(5))]
sim_model = Model(T ~ A+2B*C, Σ, LogitLink(), BinomialResponse(100))
y = simulate(sim_model, df)

# simulate two traits with two variance components using GLMM
# the two traits can have different response distribution
Σ = [VarianceComponent([0.2, 0.3], GRM),
      VarianceComponent([0.8, 0.7], eye(5))]
μ = [T1 ~ A+2B*C, T2 ~ C+log(C)+3.0]
links = [IdentityLink(), LogitLink()]
resp_dists = [NormalResponse(1.0), PoissonResponse()]
sim_model = Model(μ, Σ, links, resp_dists, )
y = simulate(sim_model, df)

# simulate two traits with two variance components with cross covariances
# using GLMM
A = [0.2 -0.1; -0.1 0.3]
B = [0.8 -0.2; -0.2 0.7]
I = eye(5)
Σ = [VarianceComponent(A, GRM),
      VarianceComponent(B, I)]
μ = [T1 ~ A+2B*C, T2 ~ C+log(C)+3.0]
sim_model = Model(μ, Σ, IdentityLink(), NormalResponse(1.0))
y = simulate(sim_model, df)

μ = [T1 ~ A+2B*C, T2 ~ C+log(C)+3.0]
sim_model = Model(μ, (@vc A ⊗ GRM + B ⊗ I), IdentityLink(),
  NormalResponse(1.0))
y = simulate(sim_model, df)

μ = Y ~ 0.2A+B+2.0
K = GRM
I = eye(5)
Σ = [VarianceComponent(0.8, K), VarianceComponent(0.2, I)]
model = Model(μ, Σ, LogLink(), PoissonResponse())
simulate(model, df)

#end
