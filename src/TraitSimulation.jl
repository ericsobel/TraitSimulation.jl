"""
This module simulates under the generalized linear model (GLM) and
generalized linear mixed model (GLMM).
"""
module TraitSimulation

export Model, simulate

using DataFrames, Distributions

"""
A type to store simulation parameters.
"""
type Model

  """
  Specify the formula of the simulation, e.g. TC ~ AGE + SNP1*SNP2 + HDL
  Using Formula of DataFrame.jl?
  """
  formula::Formula

  """
  Specify the link function, GLM.jl currently supports:
  1) CauchitLink 2) CloglogLink 3) IdentityLink 4) InverseLink
  5) LogitLink 6) LogLink 7) ProbitLink 8) SqrtLink
  """
  link::AbstractString

  """
  Specify the distribution of the response:
    1) Binomial 2) Gamma 3) Normal 4) Poisson 5) Exponential
    6) Inverse Gaussian 7) Bernoulli etc.
  """
  distribution::AbstractString

  """
  Additional parameters for the distribution, e.g. variance for normal,
  N for binomial, etc.
  """
  dist_parameters::Vector{Float64}

end

"""
Expand the right hand side of the formula
"""
function expand_rhs!(rhs::Expr, df::Symbol, df_nameset::Set{Symbol})
  for i=1:size(rhs.args,1)
    if typeof(rhs.args[i]) == Expr
      expand_rhs!(rhs.args[i], df, df_nameset)
    elseif typeof(rhs.args[i]) == Symbol && in(rhs.args[i], df_nameset)
      rhs.args[i] = parse(string(df, "[", ":", rhs.args[i], "]"))
    end
  end
end

"""
Simulate traits based on model specified in "model" using data
stored in "data_frame".
"""
function simulate(model::Model, data_frame::DataFrame)

  # parse the left and right hand side of the simulation formula
  lhs = model.formula.lhs
  rhs = model.formula.rhs

  # calculate the fixed effect component
  expand_rhs!(rhs, :data_frame, Set(names(data_frame)))
  @eval calc_rhs(data_frame)=$rhs
  fixed_eff = calc_rhs(data_frame)

end

end # end module
