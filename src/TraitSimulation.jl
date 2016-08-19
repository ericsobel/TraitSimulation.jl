"""
This module simulates under the generalized linear model (GLM) and
generalized linear mixed model (GLMM).
"""
module TraitSimulation

export Model, simulate

using DataFrame, Distributions

"""
Stores the formula of a model
"""
type ModelFormula
  lhs::Symbol
  rhs::Expr
end

"""
Construct a model formula following DataFrames convention, except that
coefficients are allowed in front of every term.
"""
macro ~(lhs, rhs)
  ex = Expr(:call, :ModelFormula, Base.Meta.quot(lhs), Base.Meta.quot(rhs))
  return ex
end

"""
A type to store simulation parameters.
"""
type Model

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
  parameters::Vector{Float64}

  """
  Specify the formula of the simulation, e.g. TC ~ AGE + SNP1*SNP2 + HDL
  Using Formula of DataFrame.jl?
  """
  formula::ModelFormula

end

"""
Simulate traits based on model specified in "model" using data
stored in "data_frame".
"""
function simulate(model::Model, data_frame::DataFrame)
  # TODO: implement this function
end

"""
An example usage of the simulate function
"""

end # end module
