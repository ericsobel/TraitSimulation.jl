# TraitSimulation.jl

A Julia module to perform trait simulation.

## Installation

TraitSimulation requires Julia 0.5. To install TraitSimulation please type
the following code in a Julia REPL:

```
Pkg.clone("https://github.com/huwenboshi/TraitSimulation.jl.git")
```

## Uninstallation

To uninstall TraitSimulation please type the following code in
a Julia REPL:

```
Pkg.rm("TraitSimulation")
```

## Key Features

1. Simulate trait under Generalized Linear Model (GLM) or Generalized Linear Mixed Model (GLMM)
2. Allow simulation of multiple correlated traits
3. Support integration with [SnpArrays](https://github.com/openmendel/SnpArrays.jl)
4. Provide convenient ways to specify the simulation model

When simulating under the GLM, traits are sampled from

$\text{Distribution}(\text{Link}^{-1}(\mu),
  \text{other distribution-specific parameters}).$
Here, \\(\mu = X\beta\\) is the mean parameter where \\(X\\) is the design matrix
and \\(\beta\\) the fixed effects, \\(\text{Link}\\) the link
function, and \\(\text{Distribution}\\) the response distribution.
Depending on the type of distribution, different distribution-specific
parameters are required. For example, simulating a normally distributed
trait with identity link requires the variance parameter \\(\sigma^2\\), i.e. one
samples from \\(N(\mu, \sigma^2)\\). To simulate a Poisson distributed trait
with log link function, one samples traits from \\(Pois(\exp(\mu))\\).

Simulating under the GLMM is similar to simulating under the GLM, except
that random effects are allowed in \\(\mu\\), e.g. \\(\mu = X\beta + Z u + \epsilon\\)
, in which \\(u\\) is the vector of random effects. Typically, \\(\mu\\) is drawn from
the normal distribution \\(N(X\beta, \sum_{i=1}^k \sigma^2_i V)\\), where
\\(\sigma^2_i\\) are the variance components. Correlated traits can be simulated
by allowing for cross covariances.

## Contact

[Huwenbo Shi](https://huwenboshi.github.io) (shihuwenbo [at] ucla [dot] edu)
