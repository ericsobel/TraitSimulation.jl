module TraitSimulationTest

using DataFrames, TraitSimulation, SnpArrays

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
sim_model = FixedEffectModel(@formula(T ~ A+2B*C),
                             IdentityLink(), NormalResponse(1.0))
y = simulate(sim_model, df)

# simulate two traits with the same link and response using GLM
μ = [@formula(T1 ~ A+2B*C), @formula(T2 ~ A+2log(B+C)+2.0)]
sim_model = FixedEffectModel(μ, IdentityLink(), NormalResponse(1.0))
y = simulate(sim_model, df)

# simulate three traits with different link and response using GLM
μ = [@formula(T1 ~ A+2B*C), @formula(T2 ~ A+2log(B+C)+2.0),
     @formula(T3 ~ A+B+C+1.0)]
links = [IdentityLink(), LogitLink(), LogLink()]
resp_dists = [NormalResponse(1.0), BinomialResponse(100), PoissonResponse()]
sim_model = FixedEffectModel(μ, links, resp_dists)
y = simulate(sim_model, df)

# simulate a single trait with two variance components using GLMM
Σ = [VarianceComponent(0.2, GRM), VarianceComponent(0.8, eye(5))]
sim_model = MixedEffectModel(@formula(T~A+2B*C), Σ, LogitLink(),
                             BinomialResponse(100))
y = simulate(sim_model, df)

# simulate two traits with two variance components using GLMM
# the two traits can have different response distribution
Σ = [VarianceComponent([0.2, 0.3], GRM),
      VarianceComponent([0.8, 0.7], eye(5))]
μ = [@formula(T1 ~ A+2B*C), @formula(T2 ~ C+log(C)+3.0)]
links = [IdentityLink(), LogitLink()]
resp_dists = [NormalResponse(1.0), PoissonResponse()]
sim_model = MixedEffectModel(μ, Σ, links, resp_dists)
y = simulate(sim_model, df)

# simulate two traits with two variance components with cross covariances
# using GLMM
A = [0.2 -0.1; -0.1 0.3]
B = [0.8 -0.2; -0.2 0.7]
I = eye(5)
Σ = [VarianceComponent(A, GRM),
     VarianceComponent(B, I)]
μ = [@formula(T1 ~ A+2B*C), @formula(T2 ~ C+log(C)+3.0)]
sim_model = MixedEffectModel(μ, Σ, IdentityLink(), NormalResponse(1.0))
y = simulate(sim_model, df)

μ = [@formula(T1 ~ A+2B*C), @formula(T2 ~ C+log(C)+3.0)]
sim_model = MixedEffectModel(μ, (@vc A ⊗ GRM + B ⊗ I), IdentityLink(),
  NormalResponse(1.0))
y = simulate(sim_model, df)

μ = @formula(Y ~ 0.2A+B+2.0)
K = GRM
I = eye(5)
Σ = [VarianceComponent(0.8, K), VarianceComponent(0.2, I)]
model = MixedEffectModel(μ, Σ, LogLink(), PoissonResponse())
simulate(model, df)

model = RandomEffectModel(:A, Σ, LogLink(), PoissonResponse())
simulate(model, df)

model = RandomEffectModel(:A, Σ, LogLink(), PoissonResponse())
simulate(model, df; pattern=0.1)

# testing SnpArrays
data = SnpArray(Pkg.dir("TraitSimulation") * "/docs/hapmap3")
sim_model = FixedEffectModel(@formula(T ~ x1+2x2*x3),
                             IdentityLink(), NormalResponse(1.0))
y = simulate(sim_model, data)

end
