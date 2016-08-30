# Model Specifications

This section introduces types for specifying the
simulation model. Other types mentioned in this section are introduced
in detail in the following sections.

## SimulationModel

```SimulationModel``` is an abstract supertype for each spefic
simulation model.

## FixedEffectModel

```FixedEffectModel``` is a type to specify simulations under the
fixed-effect model. It's a subtype of ```SimulationModel```.

```julia

FixedEffectModel(frml::FormulaType,
                 link::LinkFunctionType
                 dist::ResponseDistributionType)
```

```frml``` can be a single formula, in which case, a single trait
will be simulated, or an array of formula, in which case, multiple
traits will be simulated. Here, we use the ```Formula``` type
in the ```DataFrames``` module.

```link``` can be a single ```LinkFunction```, in which case, all traits
will be simulated under the same link function, or an array of
```LinkFunction```, in which case, each trait will be simulated under
its own link function.

```dist``` can be a single ```ResponseDistribution```, in which case,
all traits will have the same response distribution, or an array
of ```ResponseDistribution```, in which case, each trait will have
its own response distribution.

## RandomEffectModel

```RandomEffectModel``` is a type to specify simulations under the random
effect model. It's a subtype of ```SimulationModel```.

```julia
RandomEffectModel(traits::TraitType,
                  vc::Vector{VarianceComponent},
                  link::LinkFunctionType,
                  resp_dist::ResponseDistributionType)
```

```traits``` can be a single symbol, in which case, a single trait
will be simulated, or an array of symbol, in which case multiple
traits will be simulated.

```vc``` is an array of variance components / cross covariances and the
corresponding covariance matrices. See the section on the
```VarianceComponent``` type for more detail.

## MixedEffectModel

```MixedEffectModel``` is a type of specify simulations under a
fixed-effect model. It's a subtype of ```SimulationModel```.

```julia
MixedEffectModel(formula::FormulaType,
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
cross covariance matrix to be simulated. If ```var_comp``` is of type
```Vector{Float64}```, it will be interpreted as a cross covariance
matrix with off-diagnal elements equal to 0.

### operator ⊗

The TraitSimulation module implements the ```⊗``` operator to compute
the Kronecker product between two matrices.

```julia
⊗(A, B) = kron(A,B)
```

### macro @vc

The TraitSimulation module provides a macro ```@vc``` to simplify
the specification of variance components. The following code snippet
lists some equivalent ways to specify the variance components in the
simulation model. The terms in the expression to the right of ```@vc```
must be defined.

```julia
# a fake GRM
K = cor(data')
I = eye(npeople)
A = [0.2 -0.1; -0.1 0.3]
B = [0.8 -0.2; -0.2 0.7]

# the following two ways to specify Σ are equivalent
Σ = [VarianceComponent(0.2, K), VarianceComponent(0.8, I)]
Σ = @vc 0.2K + 0.8I 

# the following two ways to specify Σ are equivalent
Σ = [VarianceComponent(A, K), VarianceComponent(B, I)]
Σ = @vc A ⊗ K + B ⊗ I
```

## VarianceComponentType

```VarianceComponentType``` is a type alias to specify the
variance component or covariance matrices.

```julia
typealias VarianceComponentType
  Union{Float64, Vector{Float64}, Matrix{Float64}}
```

## FormulaType

```FormulaType``` is a type alias to specify a single formula or a
vector of formulae.

## TraitType

```TraitType``` is a type alias to specify a single symbol or a
vector of symbols.
