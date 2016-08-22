"""
This module simulates under the generalized linear model (GLM) and
generalized linear mixed model (GLMM).
"""
module TraitSimulation

export Model,
       simulate,
       NormalResponse,
       PoissonResponse,
       ExponentialResponse,
       BernoulliResponse,
       BinomialResponse,
       GammaResponse,
       InverseGaussianResponse,
       TResponse,
       CauchitLink,
       CloglogLink,
       IdentityLink,
       InverseLink,
       LogitLink,
       ProbitLink,
       SqrtLink,
       LogLink

using DataFrames,
      Distributions

"""
A list of types to store inverse link functions
TODO: Implement a "toString" function for print
"""
type CauchitLink link_inv::Function end
CauchitLink() = CauchitLink(x::Float64 -> tan(pi*(x-0.5)))

type CloglogLink link_inv::Function end
CloglogLink() = CloglogLink(x::Float64 -> 1.0-exp(-exp(x)))

type IdentityLink link_inv::Function end
IdentityLink() = IdentityLink(x::Float64 -> x)

type InverseLink link_inv::Function end
InverseLink() = InverseLink(x::Float64 -> 1.0/x)

type LogitLink link_inv::Function end
LogitLink() = LogitLink(x::Float64 -> 1.0/(1.0+exp(-x)))

const normal_dist = Normal(0.0, 1.0)
type ProbitLink link_inv::Function end
ProbitLink() = ProbitLink(x::Float64 -> pdf(normal_dist, x))

type SqrtLink link_inv::Function end
SqrtLink() = SqrtLink(x::Float64 -> x*x)

type LogLink link_inv::Function end
LogLink() = LogLink(x::Float64 -> exp(x))

const LinkFunction = Union{CauchitLink, CloglogLink,
  IdentityLink, InverseLink, LogitLink, ProbitLink,
  SqrtLink, LogLink}

"""
A list of types to store distribution parameters in simulations
TODO: Implement a "toString" function for print
"""
type NormalResponse σ::Float64 end
type PoissonResponse end
type ExponentialResponse end
type BernoulliResponse end
type BinomialResponse n::Float64 end
type GammaResponse shape::Float64 end
type InverseGaussianResponse λ::Float64 end
type TResponse ν::Float64 end
const ResponseDistribution = Union{NormalResponse, PoissonResponse,
  ExponentialResponse, BernoulliResponse, BinomialResponse,
  GammaResponse, InverseGaussianResponse, TResponse}

"""
A type to store variance component and its covariance matrix
TODO: Implement a "toString" function for print
"""
type VarianceComponent
  var_comp::Float64
  cov_mat::Matrix{Float64}
end

"""
A type to store simulation parameters.
"""
type Model

  """
  Specify the formula of the simulation, e.g. TC ~ AGE + SNP1*SNP2 + HDL
  Using Formula of DataFrame.jl?
  """
  formula::Union{Formula, Vector{Formula}}

  """
  Specify the link function, GLM.jl currently supports:
  1) CauchitLink 2) CloglogLink 3) IdentityLink 4) InverseLink
  5) LogitLink 6) LogLink 7) ProbitLink 8) SqrtLink
  """
  link::LinkFunction

  """
  Specify the distribution of the response:
    1) Binomial 2) Gamma 3) Normal 4) Poisson 5) Exponential
    6) Inverse Gaussian 7) Bernoulli etc.
  """
  resp_dist::ResponseDistribution

  """
  Specify the variance components for GLMM
  """
  vc::Vector{VarianceComponent}

  """
  Specify the cross covariance for the covariances
  """
  cross_cov::Matrix{Float64}

  """
  Specify type of distribution for variance component
  """

end

"""
Construct a Model object without random effect component
"""
function Model(formula::Union{Formula, Vector{Formula}}, link::LinkFunction,
  resp_dist::ResponseDistribution)
  return Model(formula, link, resp_dist, Vector{Float64}(), Matrix{Float64}())
end

"""
Expand the right hand side of the formula
"""
const operands = Set([:*, :/, :^])
function expand_rhs!(rhs::Expr, df::Symbol, df_nameset::Set{Symbol})
  for i=1:size(rhs.args,1)
    if typeof(rhs.args[i]) == Expr
      expand_rhs!(rhs.args[i], df, df_nameset)
    elseif typeof(rhs.args[i]) == Symbol && in(rhs.args[i], df_nameset)
      rhs.args[i] = parse(string(df, "[", ":", rhs.args[i], "]"))
    elseif typeof(rhs.args[i]) == Symbol && in(rhs.args[i], operands)
      rhs.args[i] = parse(string(".", rhs.args[i]))
    end
  end
end

"""
Simulate trait by sampling from the specified distribution
"""
function calc_trait(μ::Vector{Float64}, resp_dist::ResponseDistribution)

  if typeof(resp_dist) == NormalResponse
    return rand(Normal(0, resp_dist.σ), size(μ, 1)) + μ

  elseif typeof(resp_dist) == TResponse
    return rand(TDist(resp_dist.ν), size(μ, 1)) + μ

  elseif typeof(resp_dist) == PoissonResponse
    return map(x -> rand(Poisson(x)), μ)

  elseif typeof(resp_dist) == ExponentialResponse
    return map(x -> rand(Exponential(x)), μ)

  elseif typeof(resp_dist) == BernoulliResponse
    return map(x -> rand(Bernoulli(x)), μ)

  elseif typeof(resp_dist) == BinomialResponse
    return map(x -> rand(Binomial(resp_dist.n, x)), μ)

  elseif typeof(resp_dist) == GammaResponse
    return map(x -> rand(Gamma(resp_dist.shape, x)), μ)

  elseif typeof(resp_dist) == InverseGaussianResponse
    return map(x -> rand(InverseGaussian(x, resp_dist.λ)), μ)

  else
    # TODO: throw an exception here
    return nothing
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
  # TODO: throw an exception when formula cannot be evaluated
  expand_rhs!(rhs, :data_frame, Set(names(data_frame)))
  @eval calc_rhs(data_frame)=$rhs
  fixed_eff = calc_rhs(data_frame)

  # TODO: add random effect here
  num_indvs = size(fixed_eff, 1)
  rand_eff = zeros(Float64, num_indvs)

  # calculate the mean param: link^-1(fixed_eff + rand_eff)
  η = fixed_eff+rand_eff
  μ = map(model.link.link_inv, η)

  # simulate the trait
  y = calc_trait(μ, model.resp_dist)
  y = reshape(y, size(y,1), 1)

  # return a data frame
  y = convert(DataFrame, y)
  names!(y, [lhs])
  return y

end

end # end module
