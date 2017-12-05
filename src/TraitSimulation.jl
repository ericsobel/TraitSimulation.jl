"""
This module simulates under the generalized linear model (GLM) and
generalized linear mixed model (GLMM).
"""
module TraitSimulation

using DataFrames,
      StatsModels,
      Distributions,
      SnpArrays

export SimulationModel,
       FixedEffectModel,
       RandomEffectModel,
       MixedEffectModel,
       VarianceComponent,
       NormalResponse,
       PoissonResponse,
       ExponentialResponse,
       BernoulliResponse,
       BinomialResponse,
       GammaResponse,
       InverseGaussianResponse,
       TResponse,
       DiracResponse,
       CauchitLink,
       CloglogLink,
       IdentityLink,
       InverseLink,
       LogitLink,
       ProbitLink,
       SqrtLink,
       LogLink,
       simulate,
       @vc, ⊗

include("link_functions.jl")
include("response_distributions.jl")
include("model_specifications.jl")
include("generate_simulations.jl")

end
