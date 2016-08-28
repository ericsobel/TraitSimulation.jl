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
Define kronecker product operator
"""
⊗(A, B) = kron(A,B)

"""
Parse a variance component
"""
macro vc(expr::Expr)
  
  expr_str = string(expr)
  
  ret_expr = quote
    expr = parse($expr_str)
    [VarianceComponent(eval(expr.args[i].args[2]),
                       eval(expr.args[i].args[3]))
     for i=2:size(expr.args,1)]
  end

  return esc(ret_expr)

end

include("link_functions.jl")
include("response_distributions.jl")
include("model_specifications.jl")
include("generate_simulations.jl")

end # end module
