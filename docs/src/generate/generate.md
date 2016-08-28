# Generate Simulations

This section introduces the ```simulate``` function for simulating traits
according to the specified model.

## simulate

```julia
simulate(model::Model, data_frame::DataFrame)
```

The ```simulate``` function simulates traits according to the
```Model``` using data stored in the ```DataFrame```.

## TODO

1. add a parameter to simulate missingness
2. add other functionalities in old Mendel
3. simulate based on SnpArrays
4. simulate based on float arrays
