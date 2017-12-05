# Getting Started

This section provides some examples on how to use the
TraitSimulation module.

## Overview

The TraitSimulation module allows users to simulate phenotypes under various
models, including fixed-effect, random-effect, and mixed-effect models. The
simulation procedure takes a two-step approach. In the first step, the user
constructs a simulation model using the ```SimulationModel``` type. In the
second step, the user calls the ```simulate``` function along with the data to
generate the simulated phenotype, which is appended as a column to the input
data frame.

## Simulate a single trait

### Simulate a trait using genotype data in PLINK format

TraitSimulation module can take genotype data in PLINK format through
[SnpArrays](https://github.com/OpenMendel/SnpArrays.jl). For example, the
following code reads in genotype data from the PLINK files "hapmap3.bed",
"hapmap3.bim", "hapmap3.fam", and then simulate the trait based on the first
3 SNPs using the ```simulate```(see below for more details). More specifically,
the code below samples phenotypes from

$y \sim N(\mu, 1.0)$

where

$\mu = x_1+2x_2x_3.$


```julia
using DataFrames, Distributions, SnpArrays, TraitSimulation, StatsModels

# Load genotype data in PLINK format
data = SnpArray("hapmap3")

# Create the fixed-effect simulation model
sim_model = FixedEffectModel(@formula(T ~ x1+2x2*x3),
                             IdentityLink(), NormalResponse(1.0))

# Generate the simulation
y = simulate(sim_model, data)
```

### Simulate a trait with PC's as covariates

To include principal components as covariates in the simulation procedure,
one can apply ```pca``` in the ```SnpArrays``` module to obtain the PCs
and then include them as a column in the input data frame. The following
code simulate a trait using the first 3 SNPs and 2 PCs. More specifically,
the code below samples phenotypes from

$y \sim N(\mu, 1.0)$

where

$\mu = x_1+2x_2x_3 + PC1 + PC2.$

```julia
using DataFrames, Distributions, SnpArrays, TraitSimulation, StatsModels

# Load Genotype data in PLINK format
data = SnpArray("hapmap3")

# Compute the PC of each SNP
pcscore,_,_ = pca(data)
pcscore = convert(DataFrame, pcscore)
names!(pcscore, [Symbol("PC$i") for i in 1:size(pcscore,2)])

# Combine genotype data with PC results
data = convert(DataFrame, convert(Matrix{Float64}, data))
data = hcat(data, pcscore)

# Create the simulation model
sim_model = FixedEffectModel(@formula(T ~ x1+2x2*x3+PC1+PC2),
                             IdentityLink(), NormalResponse(1.0))

# Generate the simulations
y = simulate(sim_model, data)
```

### Simulate a trait under random effect model

The following code simulates a trait under the random effect model with
two variance components, where one of the covariance matrix is the genetic
relationship matrix (GRM) and the other the identity matrix. More specifically,
the code below samples phenotype from

$y \sim N(0, 0.8 K + 0.2I),$

where, \\(K\\) is the GRM estimated from genotype data using the ```grm```
function in ```SnpArrays```.

```julia
using DataFrames, Distributions, SnpArrays, TraitSimulation, StatsModels

# Load Genotype data in PLINK format
data = SnpArray("hapmap3")
npeople = size(data,1)

# Compute GRM using the grm function in SnpArrays
K = grm(data)
I = eye(npeople)
Σ = [VarianceComponent(0.8, K), VarianceComponent(0.2, I)]

# Create the simulation model (:T is the name of the simulated trait)
sim_model = RandomEffectModel(:T, Σ, IdentityLink(), NormalResponse(1.0))

# Generate the simulations
y = simulate(sim_model, data)
```

We also provide the ```@vc``` macro to simplify the specification of the
covariances. Note, the variables K and I must be defined before calling
the ```@vc``` macro. For example, the following two lines of code are
equivalent.

```julia
Σ = @vc 0.8K + 0.2I
Σ = [VarianceComponent(0.8, K), VarianceComponent(0.2, I)]
```

### Simulate a trait under mixed effect model

The following code simulates a trait under the mixed effect model with
two variance components as well as a fixed effect \\(\mu\\). More specifically,
the code below samples phenotype from

$y \sim N(\mu, 0.8 K + 0.2I),$

where, \\(\mu = x_1 + 2x_2x_3\\) the fixed effect vector \\(K\\) the GRM
estimated from genotype data using the ```grm``` function in ```SnpArrays```.

```julia
using DataFrames, Distributions, SnpArrays, TraitSimulation, StatsModels

# Load Genotype data in PLINK format
data = SnpArray("hapmap3")
npeople = size(data,1)

# Compute GRM using the grm function in SnpArrays
K = grm(data)
I = eye(npeople)
Σ = [VarianceComponent(0.8, K), VarianceComponent(0.2, I)]
μ = @formula(T ~ x1+2x2*x3)

# Create the simulation model (:T is the name of the simulated trait)
sim_model = MixedEffectModel(μ, Σ, IdentityLink(), NormalResponse(1.0))

# Generate the simulations
y = simulate(sim_model, data)
```

### Simulate a trait with non-normal distributions

In the previous examples, all trait values are sampled from the normal
distribution. However, users can also simulate traits with non-normal
distributions, by specifying the type of distribution and the corresponding
link function. For example, the following code simulate a trait under
the Poisson distribution. More specifically, the following code samples
phenotype from

$y \sim Pois(\lambda)$

where

$\lambda = \exp(x_1+2x_2x_3).$

```julia
using DataFrames, Distributions, SnpArrays, TraitSimulation, StatsModels

# Load genotype data in PLINK format
data = SnpArray("hapmap3")

# Create the fixed-effect simulation model
sim_model = FixedEffectModel(@formula(T ~ x1+2x2*x3),
                             LogLink(), PoissonResponse())

# Generate the simulation
y = simulate(sim_model, data)
```

## Simulate multiple traits

### Simulate three independent traits with different response distributions

The following code simulates three independent traits (\\(y_1, y_2, y_3\\))
with different response distributions.

$\mu_1 = 0.2x_1 + 3.0, y_1 \sim \text{Bin}(100, \mu_1)$

$\mu_2 = 0.1x_3 + 2.0, y_2 \sim \text{Pois}(\mu_2)$

$\mu_3 = 0.3x_4, y_3 \sim N(\mu_3, 2.0)$

```julia
using DataFrames, Distributions, SnpArrays, TraitSimulation, StatsModels

# Load genotype data in PLINK format
data = SnpArray("hapmap3")

# Create the simulation model
μ = [y1 ~ 0.2x1+3.0, y2 ~ 0.1x3+2.0, y3 ~ 0.3x4]
link = [LogitLink(), LogLink(), IdentityLink()]
dist = [BinomialResponse(100), PoissonResponse(), NormalResponse(2.0)]
sim_model = FixedEffectModel(μ, link, dist)

# Generate the simulations
y = simulate(sim_model, data)
```
### Simulate two correlated traits

Users can also simulate two correlated traits. For example, the following
code simulates two traits from

$\begin{equation}
\left(
\begin{array}{c}
y_1 \\
y_2
\end{array}
\right)
\sim
N \left(
\left(            
\begin{array}{c}
x_1 + 0.2 \\
x_3 + 0.1
\end{array}
\right),
A \otimes K + B \otimes I       
\right)
\end{equation}$

where \\(A\\) and \\(B\\) are the cross covariances.

```julia
using DataFrames, Distributions, SnpArrays, TraitSimulation, StatsModels

# Load Genotype data in PLINK format
data = SnpArray("hapmap3")
npeople = size(data,1)

# Compute GRM using the grm function in SnpArrays
K = grm(data)
I = eye(npeople)
A = [0.2 -0.1; -0.1 0.3]
B = [0.8 -0.2; -0.2 0.7]
μ = [y1 ~ x1+0.2, y2 ~ x3+0.1]

# Create the simulation model
sim_model = MixedEffectModel(μ, (@vc A ⊗ K + B ⊗ I),
                             IdentityLink(), NormalResponse(1.0))

# Generate the simulations
y = simulate(sim_model, data)
```

## Write simulations to file

```TraitSimulation.jl``` allows users to write simulated traits directly
into a file by specifying the ```out``` parameter in the ```simulate```
funciton. If the ```out``` parameter is not specified, then simulated traits
will be stored only in memory.

The following code simulates a trait and stores the output in
```output.txt``` .
```julia
using DataFrames, Distributions, SnpArrays, TraitSimulation, StatsModels

data = SnpArray("hapmap3")
data = convert(DataFrame, convert(Matrix{Float64}, data, impute=true))
sim_model = FixedEffectModel(@formula(T ~ x1+2x2*x3),
                             IdentityLink(), NormalResponse(1.0))
simulate(sim_model, data, out="output.txt")
```
