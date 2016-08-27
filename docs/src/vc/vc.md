# Variance Components

This section introduces the VarianceComponent type for specifying
the variance components in mixed model simulation.

## VarianceComponent

```julia
VarianceComponent(var_comp::Union{Float64, Vector{Float64}, Matrix{Float64}}
                  cov_mat::Matrix{Float64})
```

The ```var_comp``` parameter specifies the variance component or the
cross covariance matrix to be simulated. If ```var_comp``` is of type
```Vector{Float64}```, it will be interpreted as a cross covariance
matrix with off-diagnal elements equal to 0.

## operator ⊗

The TraitSimulation module implements the ```⊗``` operator to compute
the Kronecker product between two matrices.

```julia
⊗(A, B) = kron(A,B)
```

## macro @vc

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
