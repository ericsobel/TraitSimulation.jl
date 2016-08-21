"""
This module simulates under the generalized linear model (GLM) and
generalized linear mixed model (GLMM).
"""
module TraitSimulation

export Model, simulate, NormalResponse, IdentityLinkInv

using DataFrames, Distributions

"""
A list of types to store inverse link functions
"""
type CauchitLinkInv link_inv::Function end
CauchitLinkInv() = CauchitLinkInv(x::Float64 -> tan(pi*(x-0.5)))

type CloglogLinkInv link_inv::Function end
CloglogLinkInv() = CloglogLinkInv(x::Float64 -> 1.0-exp(-exp(x)))

type IdentityLinkInv link_inv::Function end
IdentityLinkInv() = IdentityLinkInv(x::Float64 -> x)

type InverseLinkInv link_inv::Function end
InverseLinkInv() = InverseLinkInv(x::Float64 -> 1.0/x)

type LogitLinkInv link_inv::Function end
LogitLinkInv() = LogitLinkInv(x::Float64 -> 1.0/(1.0+exp(-x)))

const normal_dist = Normal(0.0, 1.0)
type ProbitLinkInv link_inv::Function end
ProbitLinkInv() = ProbitLinkInv(x::Float64 -> pdf(normal_dist, x))

type SqrtLinkInv link_inv::Function end
SqrtLinkInv() = SqrtLinkInv(x::Float64 -> x*x)

type LogLinkInv link_inv::Function end
LogLinkInv() = LogLinkInv(x::Float64 -> exp(x))

const LinkInvFunction = Union{CauchitLinkInv, CloglogLinkInv,
  IdentityLinkInv, InverseLinkInv, LogitLinkInv, ProbitLinkInv,
  SqrtLinkInv, LogLinkInv}

"""
A list of types to store distribution parameters in simulations
"""
type NormalResponse σ::Float64 end
type PoissonResponse end
type ExponentialResponse end
type BernoulliResponse end
type BinomialResponse n::Float64 end
type GammaResponse shape::Float64 end
type InverseGaussianResponse λ::Float64 end
const ResponseDistribution = Union{NormalResponse, PoissonResponse,
  ExponentialResponse, BernoulliResponse, BinomialResponse,
  GammaResponse, InverseGaussianResponse}

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
  link::LinkInvFunction

  """
  Specify the distribution of the response:
    1) Binomial 2) Gamma 3) Normal 4) Poisson 5) Exponential
    6) Inverse Gaussian 7) Bernoulli etc.
  """
  distribution::ResponseDistribution

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
Compute the mean of the GLM or GLMM, returns a copy rather than a reference
"""

function calc_mean(η::Vector{Float64}, link::LinkInvFunction)
  return map(link.link_inv, η)
end

"""
Simulate trait by sampling from the specified distribution
"""
const supported_distributions = Set([:Normal, :Gamma, :Binomial, :Poisson,
  :Exponential, :InverseGaussian, :Bernoulli])
function calc_trait(μ::Float64, distribution::Dict)

  # check if the distribution dictionary has a name field

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

  # calculate the mean param of the distribution, i.e. calculate
  # link^-1(fixed_eff + rand_eff)
  η = fixed_eff+rand_eff
  μ = calc_mean(η, model.link)
 

end

end # end module
