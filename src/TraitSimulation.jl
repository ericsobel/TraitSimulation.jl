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
"""
# TODO: Implement a "toString" function for print
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
"""
# TODO: Implement a "toString" function for print
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
"""
# TODO: Implement a "toString" function for print
type VarianceComponent
  var_comp::Union{Float64, Vector{Float64}, Matrix{Float64}}
  cov_mat::Matrix{Float64}
end

"""
A type to store simulation parameters.
"""
# TODO: Implement a "toString" function for print
# TODO: Implement sampling from multivariate t and other elliptical
# distributions
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
  link::Union{LinkFunction, Vector{LinkFunction}}

  """
  Specify the distribution of the response:
    1) Binomial 2) Gamma 3) Normal 4) Poisson 5) Exponential
    6) Inverse Gaussian 7) Bernoulli etc.
  """
  resp_dist::Union{ResponseDistribution, Vector{ResponseDistribution}}

  """
  Specify the variance components for GLMM
  """
  vc::Vector{VarianceComponent}

end

"""
Construct a Model object without random effect component
"""
function Model(formula::Union{Formula, Vector{Formula}}, link::LinkFunction,
  resp_dist::ResponseDistribution)

  return Model(formula, link, resp_dist, Vector{VarianceComponent}())

end

"""
Expand the right hand side of the formula
"""
const operands = Set([:*, :/, :^])
function expand_rhs!(rhs::Expr, df_nameset::Set{Symbol})
  for i=1:size(rhs.args,1)
    if typeof(rhs.args[i]) == Expr
      expand_rhs!(rhs.args[i], df_nameset)

    elseif typeof(rhs.args[i]) == Symbol && in(rhs.args[i], df_nameset)
      rhs.args[i] = parse(string(:x, "[", ":", rhs.args[i], "]"))

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
Simulate random effect
"""
function calc_randeff(vc::Vector{VarianceComponent},
  npeople::Int64, ntraits::Int64)

  # create the covariance matrix
  ncomponents = length(vc)
  cov_mat = zeros(Float64, ntraits*npeople, ntraits*npeople)
  for i=1:ncomponents
    if typeof(vc[i].var_comp) == Float64 ||
       typeof(vc[i].var_comp) == Vector{Float64}
      cov_mat += kron(diag(vc[i].var_comp), vc[i].cov_mat)

    elseif typeof(vc[i].var_comp) == Vector{Float64}
      cov_mat += kron(vc[i].var_comp, vc[i].cov_mat)

    else
      #TODO: throw an exception here
      return nothing
    end
  end

  # sample the random effect
  randeff = zeros(Float64, ntraits*npeople)
  if ncomponents > 0
    dist = MvNormal(zeros(Float64, ntraits*npeople), cov_mat)
    randeff = rand(dist)
  end

  return randeff

end

"""
Simulate traits based on model specified in "model" using data
stored in "data_frame".
"""
function simulate(model::Model, data_frame::DataFrame)

  # get dimensions
  npeople = size(data_frame, 1)
  ntraits = typeof(model.formula)==Formula ? 1 : size(model.formula,1)
  formulae = typeof(model.formula)==Formula ? [model.formula] : model.formula

  # initialize traits
  y = zeros(Float64, npeople, ntraits)
  col_names = [formulae[i].lhs for i=1:ntraits]

  # evalute the formulae
  for i=1:ntraits
    lhs = formulae[i].lhs
    rhs = formulae[i].rhs

    # TODO: throw an exception when formula cannot be evaluated
    expand_rhs!(rhs, Set(names(data_frame)))
    y[:,i] = (@eval x -> $rhs)(data_frame)
  end

  # add random effect
  if length(model.vc) > 0
    randeff = calc_randeff(model.vc, npeople, ntraits)
    y += reshape(randeff, npeople, ntraits)
  end

  # apply inverse of link
  for i=1:ntraits
    y[:,i] = typeof(model.link)==Vector{Any} ?
             map(model.link[i].link_inv, y[:,i]) :
             map(model.link.link_inv, y[:,i])
  end

  # sample from the response distribution
  for i=1:ntraits
    y[:,i] = typeof(model.resp_dist)==Vector{Any} ?
             calc_trait(y[:,i], model.resp_dist[i]) :
             calc_trait(y[:,i], model.resp_dist)
  end

  # convert to data frame
  y = convert(DataFrame, y)
  names!(y, col_names)

  return y

end

end # end module
