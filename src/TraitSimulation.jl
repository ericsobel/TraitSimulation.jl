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
A list of types to store inverse link functions
"""
# TODO: Implement a "toString" function for print
abstract LinkFunction

type CauchitLink <: LinkFunction link_inv::Function end
CauchitLink() = CauchitLink(x::Float64 -> tan(pi*(x-0.5)))

type CloglogLink <: LinkFunction link_inv::Function end
CloglogLink() = CloglogLink(x::Float64 -> 1.0-exp(-exp(x)))

type IdentityLink <: LinkFunction link_inv::Function end
IdentityLink() = IdentityLink(x::Float64 -> x)

type InverseLink <: LinkFunction link_inv::Function end
InverseLink() = InverseLink(x::Float64 -> 1.0/x)

type LogitLink <: LinkFunction link_inv::Function end
LogitLink() = LogitLink(x::Float64 -> 1.0/(1.0+exp(-x)))

const normal_dist = Normal(0.0, 1.0)
type ProbitLink <: LinkFunction link_inv::Function end
ProbitLink() = ProbitLink(x::Float64 -> pdf(normal_dist, x))

type SqrtLink <: LinkFunction link_inv::Function end
SqrtLink() = SqrtLink(x::Float64 -> x*x)

type LogLink <: LinkFunction link_inv::Function end
LogLink() = LogLink(x::Float64 -> exp(x))

"""
A list of types to store distribution parameters in simulations
"""
# TODO: Implement a "toString" function for print
abstract ResponseDistribution

# normal response with standard deviation σ
type NormalResponse <: ResponseDistribution σ::Float64 end

# binomial response with n trials
type BinomialResponse <: ResponseDistribution n::Float64 end

# t response with degrees of freedom ν
type TResponse <: ResponseDistribution ν::Float64 end

# gamma response with shape parameter α
type GammaResponse <: ResponseDistribution α::Float64 end

# inverse gaussian response with shape parameter λ
type InverseGaussianResponse <: ResponseDistribution λ::Float64 end

# other distributions that do not require additional parameter
type PoissonResponse <: ResponseDistribution end
type ExponentialResponse <: ResponseDistribution end
type BernoulliResponse <: ResponseDistribution end
type DiracResponse <: ResponseDistribution end # returns itself

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
    return map(x -> rand(Gamma(resp_dist.α, x)), μ)

  elseif typeof(resp_dist) == InverseGaussianResponse
    return map(x -> rand(InverseGaussian(x, resp_dist.λ)), μ)

  elseif typeof(resp_dist) == DiracResponse
    return μ

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
      cov_mat += kron(diagm(vc[i].var_comp), vc[i].cov_mat)

    elseif typeof(vc[i].var_comp) == Matrix{Float64}
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
function simulate(model::Model, data_frame::DataFrame;
  missing::Union{Float64, Vector{Bool}}=1.0)

  # get dimensions
  npeople = size(data_frame, 1)
  ntraits = typeof(model.formula)==Formula ? 1 : size(model.formula,1)
  formulae = typeof(model.formula)==Formula ? [model.formula] : model.formula

  # initialize traits
  μ = zeros(Float64, npeople, ntraits)
  col_names = [parse(string(formulae[i].lhs)) for i=1:ntraits]

  # evalute the formulae
  for i=1:ntraits
    lhs = formulae[i].lhs
    rhs = formulae[i].rhs

    # TODO: throw an exception when formula cannot be evaluated
    expand_rhs!(rhs, Set(names(data_frame)))
    μ[:,i] = (@eval x -> $rhs)(data_frame)
  end

  # add random effect
  if length(model.vc) > 0
    randeff = calc_randeff(model.vc, npeople, ntraits)
    μ += reshape(randeff, npeople, ntraits)
  end

  # apply inverse of link
  for i=1:ntraits
    μ[:,i] = typeof(model.link)==Vector{LinkFunction} ?
             map(model.link[i].link_inv, μ[:,i]) :
             map(model.link.link_inv, μ[:,i])
  end

  # sample from the response distribution
  y = Array{Real,2}(npeople, ntraits)
  for i=1:ntraits
    y[:,i] = typeof(model.resp_dist)==Vector{ResponseDistribution} ?
             calc_trait(μ[:,i], model.resp_dist[i]) :
             calc_trait(μ[:,i], model.resp_dist)
  end

  # convert to data frame
  y = convert(DataFrame, y)
  names!(y::DataFrame, col_names)

  return y

end

end # end module

"""
Parse a variance component
"""
macro vc(expr::Expr)
  expr_str = string(expr)
  quote
    expr = parse($expr_str)
    [VarianceComponent(eval(expr.args[i].args[2]),
                       eval(expr.args[i].args[3]))
     for i=2:size(expr.args,1)]
  end
end

"""
Define kronecker product operator
"""
⊗(A, B) = kron(A,B)
