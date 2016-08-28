"""
This module simulates under the generalized linear model (GLM) and
generalized linear mixed model (GLMM).
"""
module TraitSimulation

export Model,
       simulate,
       @vc, ⊗,
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
       VarianceComponent

using DataFrames,
      Distributions



"""
Parse an expression specifying the covariance matrix of the random effects
Returns code for constructing the array of VarianceComponent
"""
macro vc(expr::Expr)
  
  # TODO: add code to check validity of expr before parseing it
  # throw an exception if the expr doesn't follow the syntax

  expr_str = string(expr)
  
  ret_expr = quote
    [VarianceComponent(eval(parse($expr_str).args[i].args[2]),
                       eval(parse($expr_str).args[i].args[3]))
     for i=2:size(parse($expr_str).args,1)]
  end

  return esc(ret_expr)

end

include("link_functions.jl")
include("response_distributions.jl")
include("model_specifications.jl")
include("generate_simulations.jl")

end # end module
