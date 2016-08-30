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

## MissingPattern

```MissingPattern``` is a type alias for some ways to indexing an array.

```julia
typealias MissingPattern
  Union{Float64, Vector{Bool}, Matrix{Bool}, BitArray{1}, BitArray{2},
        Vector{Int64}, Vector{Vector{Int64}}, UnitRange{Int64},
        Vector{UnitRange{Int64}}, StepRange{Int64, Int64},
        Vector{StepRange{Int64,Int64}}}
```
