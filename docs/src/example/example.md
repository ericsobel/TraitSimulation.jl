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

## Simulate a trait using genotype data in PLINK format

TraitSimulation module can take genotype data in PLINK format through
[SnpArrays](https://github.com/OpenMendel/SnpArrays.jl). For example, the
following code reads in genotype data from the PLINK files "hapmap3.bed",
"hapmap3.bim", "hapmap3.fam", and then simulate the trait based on the first
3 SNPs using the ```simulate```(see below for more details).

```julia
using DataFrames, Distributions, SnpArrays, TraitSimulation

# Load genotype data in PLINK format
data = SnpArray("hapmap3")

# Create the fixed-effect simulation model
sim_model = FixedEffectModel(@formula(T ~ x1+2x2*x3),
                             IdentityLink(), NormalResponse(1.0))

# Generate the simulation
y = simulate(sim_model, data)
```

## Simulate a trait with PC's as covariates

To include principal components as covariates in the simulation procedure,
one can apply ```pca``` in the ```SnpArrays``` module to obtain the PCs
and then include them as a column in the input data frame. The following
code simulate a trait using the first 3 SNPs and 2 PCs.

```julia
using DataFrames, Distributions, SnpArrays, TraitSimulation

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

## Simulate a trait under the variance component model

The following code simulates a trait with two variance components.

```julia
using DataFrames, Distributions, SnpArrays, TraitSimulation

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

## Generate random data set

The following code creates a data frame containing genotype
(x1, ..., x5) and phenotype (HDL and LDL) measurements for 10 individuals.

```julia
using DataFrames, Distributions, TraitSimulation

(people, snps) = (10, 5)
snp_data = Matrix{Float64}(people, snps)
freq = [0.2, 0.3, 0.4, 0.7, 0.5]
for i=1:snps
    snp_data[:,i] = rand(Binomial(2,freq[i]), people)
end

(hdl_data, ldl_data) = (Vector{Float64}(people), Vector{Float64}(people))
for i=1:people
    hdl_data[i] = rand(Uniform(20,80))
    ldl_data[i] = rand(Uniform(20,80))
end

data = [snp_data hdl_data ldl_data]
data = convert(DataFrame, data)

names!(data, [:x1, :x2, :x3, :x4, :x5, :HDL, :LDL])
```


## Simulate normal response

The following code simulates a trait (\\(y\\)) with normal response
(\\(\sigma = 1.0\\)) using the data frame created in [step 1](#first_step).

$\mu = -0.2x_1 + 0.1x_2 * x_5 + 0.3\log(\text{HDL})$

$y \sim N(\mu, 1.0)$

```julia
model = FixedEffectModel(y ~ -0.2x1+0.1x2*x5+0.3log(HDL), IdentityLink(), NormalResponse(1.0))
simulate(model, data)
```

## Simulate multiple traits

The following code simulates three traits (\\(y_1, y_2, y_3\\)) with
normal response (\\(\sigma = 1.0\\)) but different means, using the data
frame created in [step 1](#first_step).

```julia
model = FixedEffectModel([y1 ~ 3.0+0.2x1, y2 ~ 2.0+0.1x3, y3 ~ 0.3x4+HDL], IdentityLink(), NormalResponse(1.0))
simulate(model, data)
```

The following code simulates three traits (\\(y_1, y_2, y_3\\)) with different
response distributions using the data frame created in [step 1](#first_step).

$\mu_1 = 0.2x_1 + 3.0, y_1 \sim \text{Bin}(100, \mu_1)$

$\mu_2 = 0.1x_3 + 2.0, y_2 \sim \text{Pois}(\mu_2)$

$\mu_3 = 0.3x_4 + \text{HDL}, y_3 \sim N(\mu_3, 2.0)$

```julia
μ = [y1 ~ 0.2x1+3.0, y2 ~ 0.1x3+2.0, y3 ~ 0.3x4+HDL]
link = [LogitLink(), LogLink(), IdentityLink()]
dist = [BinomialResponse(100), PoissonResponse(), NormalResponse(2.0)]
model = FixedEffectModel(μ, link, dist)
simulate(model, data)
```
## Simulate random effects

The following code simulates a trait with Poisson response with
two variance components, using the data frame created in [step 1](#first_step).

$\mu = (0.2x_1 + 2.0) + x u + \epsilon, u \sim N(0, 0.04K), \epsilon \sim N(0, 0.8I), y \sim \text{Pois}(\mu)$

```julia
# a fake GRM
K = cor(snp_data')
I = eye(people)
Σ = [VarianceComponent(0.2, K), VarianceComponent(0.8, I)]
μ = y ~ 0.2x1+2.0
model = MixedEffectModel(μ, Σ, LogLink(), PoissonResponse())
simulate(model, data)
```

We also provide the ```@vc``` macro to simplify the specification of the covariances.
Note, the variables K and I must be defined before calling the ```@vc``` macro.

```julia
Σ = @vc 0.2K + 0.8I
model = MixedEffectModel(μ, Σ, LogLink(), PoissonResponse())
simulate(model, data)
```

The following code snippet simulates two traits where the random effects have cross covariances.

```julia
# a fake GRM
K = cor(data')
I = eye(npeople)
A = [0.2 -0.1; -0.1 0.3]
B = [0.8 -0.2; -0.2 0.7]
μ = [y1 ~ x1+0.2, y2 ~ x3+0.1log(HDL)+0.1]
model = MixedEffectModel(μ, (@vc A ⊗ K + B ⊗ I), IdentityLink(), NormalResponse(1.0))
simulate(model, data)
```
