# Simulating the Model

This section introduces the ```Model``` type for specifying the
simulation model and the ```simulate``` function for generating simulations
according to the model. Types mentioned in this section are introduced
in detail in the following sections.

## Model

The ```Model``` type provides two types of constructors for specifying
a simulation module with / without random effect

### Specify a simulation model with fixed-effect only

```julia

Model(frml::Union{Formula, Vector{Formula}},
      link::Union{LinkFunction, Vector{LinkFunction}},
      dist::Union{ResponseDistribution, Vector{ResponseDistribution}})
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

### Specifying the variance components

```julia
Model(frml::Union{Formula, Vector{Formula}},
      vc::Vector{VarianceComponent},
      link::Union{LinkFunction, Vector{LinkFunction}},
      dist::Union{ResponseDistribution, Vector{ResponseDistribution}})
```

In addition to the three parameters in the fixed-effect model constructor,
the mixed-effect model constructors requires the additional parameter
```vc```, an array of variance components.

## simulate

```julia
simulate(model::Model, data_frame::DataFrame)
```

The ```simulate``` function simulates traits according to the
```Model``` using data stored in the ```DataFrame```.
