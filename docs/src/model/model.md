# Model Specifications

This section introduces types for specifying the
simulation model. Other types mentioned in this section are introduced
in detail in the following sections.

## SimulationModel

```SimulationModel``` is an abstract supertype for each specific
simulation model.

## FixedEffectModel

```FixedEffectModel``` specifies simulations under the
fixed-effect model. It's a subtype of ```SimulationModel```.

```julia

FixedEffectModel(frml::FormulaLike,
                 link::LinkFunctionType
                 dist::ResponseDistributionType)
```

```frml``` can be a single formula, in which case, a single trait
will be simulated, or an array of formulae, in which case, multiple
traits will be simulated. Here, we use the ```Formula``` type
in the ```DataFrames``` module.

```link``` can be a single link function, in which case, all traits
will be simulated under the same link function, or a vector of
link functions, in which case, each trait will be simulated under
its own link function.

```dist``` can be a single response distribution, in which case,
all traits will have the same response distribution, or a vector
of response distributions, in which case, each trait will have
its own response distribution.

## RandomEffectModel

```RandomEffectModel``` specifies simulations under the random
effect model. It's a subtype of ```SimulationModel```.

```julia
RandomEffectModel(traits::TraitType,
                  vc::Vector{VarianceComponent},
                  link::LinkFunctionType,
                  resp_dist::ResponseDistributionType)
```

```traits``` can be a single symbol, in which case, a single trait
will be simulated, or a vector of symbols, in which case multiple
traits will be simulated.

```vc``` is an array of variance components / cross covariances and the
corresponding covariance matrices. See the section on the
```VarianceComponent``` type for more detail.

## MixedEffectModel

```MixedEffectModel``` is a type of specify simulations under a
fixed-effect model. It's a subtype of ```SimulationModel```.

```julia
MixedEffectModel(formula::FormulaLike,
                 vc::Vector{VarianceComponent},
                 link::LinkFunctionType,
                 resp_dist::ResponseDistributionType)
```

## VarianceComponent

```julia
VarianceComponent(var_comp::Union{Float64, Vector{Float64}, Matrix{Float64}}
                  cov_mat::Matrix{Float64})
```

The ```var_comp``` parameter specifies the variance component or the
covariance matrix to be simulated. If ```var_comp``` is of type
```Vector{Float64}```, it will be interpreted as a diagonal covariance
matrix with off-diagnal elements equal to 0.

### operator ⊗

The TraitSimulation module implements the ```⊗``` operator to compute
the Kronecker product between two matrices.

```julia
⊗(A, B) = kron(A,B)
```

### macro @vc

The TraitSimulation module provides a macro ```@vc``` to simplify
the specification of variance components. The following code
lists some examples to specify the variance components in a
simulation model. The terms in the expression to the right of ```@vc```
must be defined.

```julia
# a fake GRM
K = cor(snp_data')
I = eye(people)
A = [0.2 -0.1; -0.1 0.3]
B = [0.8 -0.2; -0.2 0.7]

Σ = @vc 0.2K + 0.8I
Σ = @vc A ⊗ K + B ⊗ I
```

## VarianceComponentType

```VarianceComponentType``` is a type alias to specify the
variance component or covariance matrices.

```julia
const VarianceComponentType =
  Union{Float64, Vector{Float64}, Matrix{Float64}}
```

## FormulaLike

```FormulaLike``` is a type alias to specify a single formula or a
vector of formulae.

## TraitType

```TraitType``` is a type alias to specify a single symbol or a
vector of symbols.
