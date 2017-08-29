# Generate Simulations

This section introduces the ```simulate``` function for simulating traits
according to the specified model.

## simulate

```julia
simulate(model::Model, data_frame::DataFrame; pattern::MissingPattern=1.0)
```

The ```simulate``` function simulates traits according to the
```Model``` using data stored in the ```DataFrame```. One can use
the ```pattern``` keyword to specify desired missing pattern. See the 
section ```MissingPattern``` for details on supported missing patterns.

## InputDataType

```InputDataType``` is a type alias for different types of input data, which
can be either a ```DataFrame``` or ```SnpArray{2}```.

## MissingPattern

```MissingPattern``` is a type alias for different ways to specify missing
entries of simulated trait(s). Entries marked as missing have the
value ```NA```.

```julia
typealias MissingPattern
  Union{Float64, Vector{Bool}, Matrix{Bool}, BitArray{1}, BitArray{2},
        Vector{Int64}, Vector{Vector{Int64}}, UnitRange{Int64},
        Vector{UnitRange{Int64}}, StepRange{Int64, Int64},
        Vector{StepRange{Int64,Int64}}}
```

### Float64

Missing pattern specified as float will be interpreted as missing rate.
For simulations involving multiple traits, the same missing rate is applied
across all traits.

### Vector{Bool}, BitArray{1}, or Vector{Int64}

Missing pattern specified as ```Vector{Bool}```, ```BitArray{1}```
or ```Vector{Int64}``` will be interpreted as an index set containing indices
of the missing entries. For simulations involving multiple traits, the same
missing pattern is applied across all traits.

### Matrix{Bool}, BitArray{2}, or Vector{Vector{Int64}}

Missing pattern specified as ```Matrix{Bool}``` or ```BitArray{2}```  with
dimension $people \times traits$ or ```Vector{Vector{Int64}}``` with length
$traits$ will be interpreted as index sets containing indices of the missing
entries. This allows for different traits to have different index sets for
missing entries.

### UnitRange{Int64} or StepRange{Int64, Int64}

Missing pattern spcified as ```UnitRange{Int64}``` or
```StepRange{Int64, Int64}``` can be used to specify a range of missing
entries. For simulations involving multiple traits, the same missing pattern
will be applied across all traits.

### Vector{UnitRange{Int64}} or Vector{StepRange{Int64, Int64}}

Missing pattern specified as ```Vector{UnitRange{Int64}}``` or
```Vector{StepRange{Int64, Int64}}``` can be used to specify different
ranges of missing entries for different traits for simulations involing
multiple traits.

