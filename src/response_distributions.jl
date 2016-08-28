"""
A list of types to store distribution parameters in simulations
"""

# super type of all response distribution types
abstract ResponseDistribution

# normal response with standard deviation σ
type NormalResponse <: ResponseDistribution
  σ::Float64
end

# binomial response with n trials
type BinomialResponse <: ResponseDistribution
  n::Float64
end

# t response with degrees of freedom ν
type TResponse <: ResponseDistribution
  ν::Float64
end

# gamma response with shape parameter α
type GammaResponse <: ResponseDistribution
  α::Float64
end

# inverse gaussian response with shape parameter λ
type InverseGaussianResponse <: ResponseDistribution
  λ::Float64
end

# Poisson distributed response
type PoissonResponse <: ResponseDistribution
end

# exponential distributed response
type ExponentialResponse <: ResponseDistribution
end

# Bernoulli distributed response
type BernoulliResponse <: ResponseDistribution
end

# Dirac distributed response
type DiracResponse <: ResponseDistribution
end
