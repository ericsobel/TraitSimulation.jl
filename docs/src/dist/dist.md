# Response Distribution

This sections lists the response distributions supported by the TraitSimulation module.

## ResponseDistribution

```ResponseDistribution``` is an abstract type. It's the super type
of each specific response distribution listed below. All sub
types of ```ResponseDistribution``` have constructors of the form
```[distribution type](...)```. Parameter(s) passed to the constructor
of each ResponseDistribution are specific to the type of the
distribution.

## NormalResponse

Constructs a normal response type, through ```NormalResponse(σ)```,
where ```σ``` is the standard deviation of the normal distribution.

## BinomialResponse

Constructs a binomial response type, through ```BinomialResponse(n)```,
where ```n``` is the number of trials in the binomial distribution.

## TResponse

Constructs a T response type, through ```TResponse(ν)```, 
where ```ν``` is the degrees of freedom. The simulated trait is a
shifted (by $\mu$) version of T response.

## GammaResponse

Constructs a Gamma response type, through ```GammaResponse(α)```,
where ```α``` is the shape parameter.

## InverseGaussianResponse

Constructs an inverse Gaussian response, through
```InverseGaussianResponse(λ)``` where ```λ``` is the shape parameter.

## PoissonResponse

Constructs a Poisson response, through ```Poisson()```.

## ExponentialResponse

Constructs an exponential response, through ```ExponentialResponse()```.

## BernoulliResponse

Constructs a Bernoulli response, through ```Bernoulli()```.

## DirectResponse

Directly output the mean, without imposing any distribution.
