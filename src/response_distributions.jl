"""
A list of types to store distribution parameters in simulations
"""


# super type of all response distribution types
abstract ResponseDistribution


# normal response with standard deviation σ
type NormalResponse <: ResponseDistribution
  σ::Float64

  # to make sure σ is greater than 0
  function NormalResponse(σ::Float64)
    if(σ <= 0)
      error("σ must be greater than zero in a normal distribution!")
    else
      new(σ)
    end
  end

end


# binomial response with n trials
type BinomialResponse <: ResponseDistribution
  n::Int64
end


# t response with degrees of freedom ν
type TResponse <: ResponseDistribution
  ν::Float64
end


# gamma response with shape parameter α
type GammaResponse <: ResponseDistribution
  α::Float64

  # to make sure α is greater than 0
  function GammaResponse(α::Float64)
    if(α <= 0)
      error("α must be greater than zero in a Gamma distribution!")
    else
      new(α)
    end
  end

end


# inverse gaussian response with shape parameter λ
type InverseGaussianResponse <: ResponseDistribution
  λ::Float64

  # to make sure λ is greater than 0
  function InverseGaussianResponse(α::Float64)
    if(λ <= 0)
      error("λ must be greater than zero in a inverse Gaussian distribution!")
    else
      new(α)
    end
  end

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


# Direct response, simply output the mean
type DirectResponse <: ResponseDistribution
end
