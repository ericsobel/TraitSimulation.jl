"""
This file contains functions to simulate traits.
"""

# define some type alias for clearer code
const MissingPattern =
  Union{Float64, Vector{Bool}, Matrix{Bool}, BitArray{1}, BitArray{2},
        Vector{Int64}, Vector{Vector{Int64}}, UnitRange{Int64},
        Vector{UnitRange{Int64}}, StepRange{Int64, Int64},
        Vector{StepRange{Int64,Int64}}}

const InputDataType = Union{DataFrame, SnpArray{2}}


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
      rhs.args[i] = parse(string(:input_data_from_user_qJvsFLOpUt,
        "[", ":", rhs.args[i], "]"))

    elseif typeof(rhs.args[i]) == Symbol && in(rhs.args[i], operands)
      rhs.args[i] = parse(string(".", rhs.args[i]))
    end
  end
end


"""
Simulate trait by sampling from the specified distribution
"""
function calc_trait(μ::Vector{Float64}, resp_dist::ResponseDistribution,
  nvarcomponents::Int64)

  if typeof(resp_dist) == NormalResponse
    if (nvarcomponents == 0)
      return rand(Normal(0, resp_dist.σ), size(μ, 1)) + μ
    else
      return μ
    end

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

  # get number of components
  nvarcomponents = length(vc)

  # normal distribution with mean 0 and variance 1
  dist = Normal(0.0, 1.0)
  randeff = zeros(Float64, npeople*ntraits)

  # construct the random effect iteratively
  for i=1:nvarcomponents

    # in case of float or vector, make them diagonal
    if typeof(vc[i].var_comp) == Float64 ||
       typeof(vc[i].var_comp) == Vector{Float64}

      cross_cov = diagm(vc[i].var_comp)
      cov_mat = vc[i].cov_mat

    # input is a matrix, no change
    else typeof(vc[i].var_comp) == Matrix{Float64}

      cross_cov = vc[i].var_comp
      cov_mat = vc[i].cov_mat

    end

    # cholesky of the kronecker product
    L = chol(cross_cov)' ⊗ chol(cov_mat)'

    # increment the cholesky
    tmp = rand(dist, npeople*ntraits)
    randeff += L * tmp

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
      missing_idx = randperm(npeople)[1:num_missing]
      df[missing_idx,i] = missing
    end

  # if missing pattern is a bool vector, bit array, step range
  elseif typeof(pattern)==Vector{Bool} || typeof(pattern)==BitArray{1} ||
         typeof(pattern)==Vector{Int64} || typeof(pattern)==UnitRange{Int64} ||
         typeof(pattern) == StepRange{Int64,Int64}
    for i=1:ntraits
      df[pattern,i] = missing
    end

  # if missing pattern is a matrix of bool or bit array
  elseif typeof(pattern) == Matrix{Bool} || typeof(pattern) == BitArray{2}
    df[pattern] = missing

  # if missing pattern is a vector of vector of int 64 or vector of ranges
  else
    for i=1:ntraits
      df[pattern[i],i] = missing
    end
  end

end


"""
Simulate traits based on model specified in "model" using data
stored in "data".
"""
function simulate(model::SimulationModel, data::InputDataType;
  pattern::MissingPattern=0.0, out::String="")

  # if data is SnpArray, convert to data frame first
  if typeof(data) == SnpArray{2}
    data = convert(DataFrame, convert(Matrix{Float64}, data, impute=true))
  end

  global input_data_from_user_qJvsFLOpUt = data

  # get dimensions
  npeople, ntraits = (size(data, 1), size(model))

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
      rhs = formulae[i].rhs
      try
        expand_rhs!(rhs, Set(names(data)))
        μ[:,i] = eval(rhs)
      catch
        throw(ErrorException(string("Failed to evaluate: ", rhs)))
      end
    end
  end

  # add random effect
  if typeof(model) == RandomEffectModel || typeof(model) == MixedEffectModel
    randeff = calc_randeff(model.vc, npeople, ntraits)
    μ += reshape(randeff, npeople, ntraits)
    nvarcomponents = length(model.vc)
  else
    nvarcomponents = 0
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
             calc_trait(μ[:,i], model.resp_dist[i], nvarcomponents) :
             calc_trait(μ[:,i], model.resp_dist, nvarcomponents)
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

  if out != ""
    CSV.write(out, y; delim='\t')
  end

  return y

end
