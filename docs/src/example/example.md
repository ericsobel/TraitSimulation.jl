# Getting Started

This section provides some examples on how to use the
TraitSimulation module.

## <a name="first_step">Random test data set</a> 

The following code creates a data frame containing genotype
(X1, ..., X5) and phenotype (HDL and LDL) measurements for 10 individuals.

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
data_frame = convert(DataFrame, data)
names!(data_frame, [:X1, :X2, :X3, :X4, :X5, :HDL, :LDL])
```
## Simulate normal response

The following code simulates a trait ($Y$) with normal response
(\\(\sigma = 1.0\\)) using the data frame created in [step 1](#first_step).

$\mu = -0.2X_1 + 0.1X_2 * X_5 + 0.3\log(\text{HDL})$

$Y \sim N(\mu, 1.0)$

```julia
model = FixedEffectModel(Y ~ -0.2X1+0.1X2*X5+0.3log(HDL), IdentityLink(), NormalResponse(1.0))
simulate(model, data_frame)
```

## Simulate multiple traits

The following code simulates three traits ($Y_1, Y_2, Y_3$) with
normal response ($\sigma = 1.0$) but different means, using the data
frame created in [step 1](#first_step).

```julia
model = FixedEffectModel([Y1 ~ 3.0+0.2X1, Y2 ~ 2.0+0.1X3, Y3 ~ 0.3X4+HDL], IdentityLink(), NormalResponse(1.0))
simulate(model, data_frame)
```

The following code simulates three traits ($Y_1, Y_2, Y_3$) with different
response distributions using the data frame created in [step 1](#first_step).

$\mu_1 = 0.2X_1 + 3.0, Y_1 \sim \text{Bin}(100, \mu_1)$

$\mu_2 = 0.1X_3 + 2.0, Y_2 \sim \text{Pois}(\mu_2)$

$\mu_3 = 0.3X_4 + \text{HDL}, Y_3 \sim N(\mu_3, 2.0)$

```julia
μ = [Y1 ~ 0.2X1+3.0, Y2 ~ 0.1X3+2.0, Y3 ~ 0.3X4+HDL]
link = [LogitLink(), LogLink(), IdentityLink()]
dist = [BinomialResponse(100), PoissonResponse(), NormalResponse(2.0)]
model = FixedEffectModel(μ, link, dist)
simulate(model, data_frame)
```
## Simulate random effects

The following code simulates a trait with Poisson response with
two variance components, using the data frame created in [step 1](#first_step).

$\mu = (0.2X_1 + 2.0) + X u + \epsilon, u \sim N(0, 0.04K), \epsilon \sim N(0, 0.8I), Y \sim \text{Pois}(\mu)$

```julia
# a fake GRM
K = cor(snp_data')
I = eye(people)
Σ = [VarianceComponent(0.2, K), VarianceComponent(0.8, I)]
μ = Y ~ 0.2X1+2.0
model = MixedEffectModel(μ, Σ, LogLink(), PoissonResponse())
simulate(model, data_frame)
```

We also provide the ```@vc``` macro to simplify the specification of the covariances.
Note, the variables K and I must be defined before calling the ```@vc``` macro.

```julia
Σ = @vc 0.2K + 0.8I
model = MixedEffectModel(μ, Σ, LogLink(), PoissonResponse())
simulate(model, data_frame)
```

The following code snippet simulates two traits where the random effects have cross covariances.

```julia
# a fake GRM
K = cor(data')
I = eye(npeople)
A = [0.2 -0.1; -0.1 0.3]
B = [0.8 -0.2; -0.2 0.7]
μ = [Y1 ~ X1+0.2, Y2 ~ X3+0.1log(HDL)+0.1]
model = MixedEffectModel(μ, (@vc A ⊗ K + B ⊗ I), IdentityLink(), NormalResponse(1.0))
simulate(model, data_frame)
```
