"""
A type to store variance component and its covariance matrix
"""
# TODO: Implement a "toString" function for print
type VarianceComponent

  """
  Stores a single variance component, a vector of variance component
  for a number of traits, or a cross covariance matrix
  """
  var_comp::Union{Float64, Vector{Float64}, Matrix{Float64}}

  """
  Stores the covariance matrix
  """
  cov_mat::Matrix{Float64}

end

"""
A type to store simulation parameters.
"""
# TODO: Implement a "toString" function for print
type Model

  """
  Specify the formula of the simulation, using DataFrames' Formula
  e.g. TC ~ AGE + 1.5SNP1*SNP2 + 2.0HDL
  """
  formula::Union{Formula, Vector{Formula}}

  """
  Specify the variance components for GLMM
  """
  vc::Vector{VarianceComponent}
  
  """
  Specify the link function, GLM.jl currently supports:
  1) CauchitLink 2) CloglogLink 3) IdentityLink 4) InverseLink
  5) LogitLink 6) LogLink 7) ProbitLink 8) SqrtLink
  """
  link::Union{LinkFunction, Vector{LinkFunction}}

  """
  Specify the distribution of the response:
    1) Binomial 2) Gamma 3) Normal 4) Poisson 5) Exponential
    6) Inverse Gaussian 7) Bernoulli etc.
  """
  resp_dist::Union{ResponseDistribution, Vector{ResponseDistribution}}

end

"""
Construct a Model object without random effect component
"""
Model(formula::Union{Formula, Vector{Formula}},
  link::Union{LinkFunction, Vector{LinkFunction}},
  resp_dist::Union{ResponseDistribution, Vector{ResponseDistribution}}) =
Model(formula, Vector{VarianceComponent}(), link, resp_dist)

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
