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
  link::Symbol

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
  dist_params::Vector{Float64}

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
const normal_dist = Normal(0.0, 1.0)
function calc_mean(η::Vector{Float64}, link::Symbol)
  if link == :IdentityLink
    return η
  elseif link == :LogitLink
    @eval logit_inv(x::Float64)=1.0/(1.0+exp(-x))
    return map(logit_inv::Float64, η)
  elseif link == :LogLink
    return map(exp, η)
  elseif link == :InverseLink
    @eval inverse_inv(x::Float64)=1.0/x
    return map(inverse_inv, η)
  elseif link == :SqrtLink
    @eval sqrt_inv(x::Float64)=x*x
    return map(sqrt_inv, η)
  elseif link == :ProbitLink
    @eval probit_inv(x::Float64)=pdf(normal_dist, x)
    return map(probit_inv, η)
  elseif link == :CauchitLink
    @eval cauchit_inv(x::Float64)=tan(pi*(x-0.5))
    return map(cauchit_inv, η)
  elseif link == :CloglogLink
    @eval cloglog_inv(x::Float64)=1.0-exp(-exp(x))
    return map(cloglog_inv, η)
  else
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
