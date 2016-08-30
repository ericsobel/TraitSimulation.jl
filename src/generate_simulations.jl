"""
This file contains functions to simulate traits.
"""


# define some some type alias for clearer code
typealias MissingPattern
  Union{Float64, Vector{Bool}, Matrix{Bool}, BitArray{1}, BitArray{2},
        Vector{Int64}, Vector{Vector{Int64}}, UnitRange{Int64},
        Vector{UnitRange{Int64}}, StepRange{Int64, Int64},
        Vector{StepRange{Int64,Int64}}}


"""
Define kronecker product operator for the symbol ⊗
"""
⊗(A, B) = kron(A,B)


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

  elseif typeof(resp_dist) == DirectResponse
    return μ

  else
    # throw an exception if the response distribution is not supported
    throw(ErrorException("Response distribution not supported."))
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
    
    # if a float or vector of float, make a diagonal matrix
    if typeof(vc[i].var_comp) == Float64 ||
       typeof(vc[i].var_comp) == Vector{Float64}
      cov_mat += kron(diagm(vc[i].var_comp), vc[i].cov_mat)

    else typeof(vc[i].var_comp) == Matrix{Float64}
      cov_mat += kron(vc[i].var_comp, vc[i].cov_mat)
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
Create some missingness in the simulated data
"""
function missing!(df::DataFrame, pattern::MissingPattern)
  
  # get dimensions
  ntraits = size(df, 2)
  npeople = size(df, 1)

  # if missing pattern is a float, sample missing entries
  if typeof(pattern) == Float64

    # sample some missing entries
    num_missing = convert(Int64, floor(pattern*npeople))
    
    for i=1:ntraits
      missing_idx = randperm(npeople)[1:npeople]
      df[missing_idx,i] = NA
    end

  # if missing pattern is a bool vector, bit array, step range
  elseif typeof(pattern)==Vector{Bool} || typeof(pattern)==BitArray{1} ||
         typeof(pattern)==Vector{Int64} || typeof(pattern)==UnitRange{Int64} ||
         typeof(pattern) == StepRange{Int64,Int64}
    for i=1:ntraits
      df[pattern,i] = NA
    end
  
  # if missing pattern is a matrix of bool or bit array
  elseif typeof(pattern) == Matrix{Bool} || typeof(pattern) == BitArray{2}
    df[pattern] = NA
  
  # if missing pattern is a vector of vector of int 64 or vector of ranges
  else
    for i=1:ntraits
      df[pattern[i],i] = NA
    end
  end
   
end


"""
Simulate traits based on model specified in "model" using data
stored in "data_frame".
"""
function simulate(model::SimulationModel, data_frame::DataFrame;
  pattern::MissingPattern=0.0)

  # get dimensions
  npeople = size(data_frame, 1)
  ntraits = size(model)

  # for simulating under fixed and mixed effect model only
  if typeof(model) == FixedEffectModel || typeof(model) == MixedEffectModel
    formulae = typeof(model.formula)==Formula?[model.formula]:model.formula
  end

  # initialize traits
  μ = zeros(Float64, npeople, ntraits)
  if typeof(model) == FixedEffectModel || typeof(model) == MixedEffectModel
    col_names = [parse(string(formulae[i].lhs)) for i=1:ntraits]
  else
    col_names = typeof(model.traits)==Symbol ? [model.traits] : model.traits
  end

  # evalute the formulae for fixed and mixed effect model
  if typeof(model) == FixedEffectModel || typeof(model) == MixedEffectModel
    for i=1:ntraits
      lhs = formulae[i].lhs
      rhs = formulae[i].rhs

      # eval the rhs, throw an exception when formula cannot be evaluated
      expand_rhs!(rhs, Set(names(data_frame)))
      try
        μ[:,i] = (@eval x -> $rhs)(data_frame)
      catch
        throw(ErrorException(string("Failed to evaluate: ", rhs)))
      end

    end
  end

  # add random effect
  if typeof(model) == RandomEffectModel || typeof(model) == MixedEffectModel
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
  for i=1:ntraits
    y[:,i] = convert(DataArray, y[:,i])
  end

  # impose missingness
  try
    missing!(y, pattern)
  catch
    throw(ErrorException("An error occured while imposing missingness."))
  end

  return y

end
