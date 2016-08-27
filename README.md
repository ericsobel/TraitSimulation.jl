##TraitSimulation

*TraitSimulation* is a julia package that provides utilities to simulate
phenotypes under a generalized linear model or generalized linear mixed model.

[![Build Status](https://travis-ci.org/huwenboshi/TraitSimulation.jl.svg?branch=master)](https://travis-ci.org/huwenboshi/TraitSimulation.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://huwenboshi.github.io/TraitSimulation.jl)

##Installation

To install *TraitSimulation*, please type the following command in a julia
interactive terminal
```julia
Pkg.clone("https://github.com/huwenboshi/TraitSimulation.jl.git")
```

##Examples


#### First, simulate some data


```julia
include("../src/TraitSimulation.jl");
using DataFrames, Distributions, TraitSimulation;
srand(1);
npeople, nsnp = (10, 5);
snp_data = Matrix{Float64}(npeople, nsnp);
freqs = [0.2, 0.3, 0.4, 0.7, 0.5];
for i=1:nsnp
    snp_data[:,i] = rand(Binomial(2,freqs[i]), npeople);
end
hdl_data, ldl_data = (Vector{Float64}(npeople),
    Vector{Float64}(npeople));
for i=1:npeople
    hdl_data[i] = rand(Uniform(20,80));
    ldl_data[i] = rand(Uniform(20,80));
end
data = [snp_data hdl_data ldl_data]
data_frame = convert(DataFrame, data);
names!(data_frame, [:X1, :X2, :X3, :X4, :X5, :HDL, :LDL])
```




<table class="data-frame"><tr><th></th><th>X1</th><th>X2</th><th>X3</th><th>X4</th><th>X5</th><th>HDL</th><th>LDL</th></tr><tr><th>1</th><td>0.0</td><td>1.0</td><td>1.0</td><td>2.0</td><td>1.0</td><td>34.162000739722814</td><td>40.79102085151763</td></tr><tr><th>2</th><td>1.0</td><td>0.0</td><td>0.0</td><td>1.0</td><td>1.0</td><td>38.76241810016405</td><td>20.474557003433645</td></tr><tr><th>3</th><td>0.0</td><td>0.0</td><td>2.0</td><td>2.0</td><td>0.0</td><td>49.316769804770075</td><td>32.658092129512156</td></tr><tr><th>4</th><td>0.0</td><td>0.0</td><td>2.0</td><td>1.0</td><td>1.0</td><td>77.11498039014404</td><td>79.99427953391682</td></tr><tr><th>5</th><td>2.0</td><td>0.0</td><td>2.0</td><td>1.0</td><td>1.0</td><td>35.09973098191831</td><td>79.19998201392798</td></tr><tr><th>6</th><td>0.0</td><td>1.0</td><td>0.0</td><td>0.0</td><td>2.0</td><td>53.34506523947434</td><td>46.22647847657751</td></tr><tr><th>7</th><td>0.0</td><td>0.0</td><td>1.0</td><td>2.0</td><td>1.0</td><td>45.48307102970789</td><td>66.39338290744263</td></tr><tr><th>8</th><td>2.0</td><td>1.0</td><td>0.0</td><td>1.0</td><td>0.0</td><td>36.871413937143785</td><td>32.568342391884244</td></tr><tr><th>9</th><td>0.0</td><td>1.0</td><td>0.0</td><td>1.0</td><td>2.0</td><td>35.0827525875335</td><td>21.222492122760034</td></tr><tr><th>10</th><td>2.0</td><td>0.0</td><td>1.0</td><td>0.0</td><td>1.0</td><td>37.26209073654137</td><td>71.57072816525965</td></tr></table>



#### Simulate a single trait with Normal response

$\mu = -0.2X_1 + 0.1X_2 \times X_5 + 0.3\log(\text{HDL} + \text{LDL})$

$y \sim N(\mu, 1.0)$


```julia
model = Model(Y ~ -0.2X1+0.1X2*X5+0.3log(HDL+LDL),
    IdentityLink(), NormalResponse(1.0))
simulate(model, data_frame)
```

#### Simulate three traits with different mean but same response distribution

$\mu_1 = 0.2X1+3.0$, $\mu_2 = 0.3X_3+2.0$, $\mu_3 = 0.3X_4+\text{HDL}$

$y_1 \sim N(\mu_1, 1.0)$, $y_2 \sim N(\mu_2, 1.0)$, $y_3 \sim N(\mu_3, 1.0)$


```julia
model = Model([Y1 ~ 0.2X1+3.0, Y2 ~ 0.1X3+2.0, Y3 ~ 0.3X4+HDL],
    IdentityLink(), NormalResponse(1.0))
simulate(model, data_frame)
```

#### Simulate three traits with Binomial, Poisson, and Normal response

$\mu_1 = 0.2X_1 + 3.0$, 
$y_1 \sim \text{Bin}(100, \mu_1)$

$\mu_2 = 0.1X_3 + 2.0$, 
$y_2 \sim \text{Pois}(\mu_2)$

$\mu_3 = 0.3X_4 + HDL$, 
$y_3 \sim N(\mu_3, 2.0)$


```julia
μ = [Y1 ~ 0.2X1+3.0, Y2 ~ 0.1X3+2.0, Y3 ~ 0.3X4+HDL]
link = [LogitLink(), LogLink(), IdentityLink()]
dist = [BinomialResponse(100), PoissonResponse(), NormalResponse(2.0)]
model = Model(μ, link, dist)
simulate(model, data_frame)
```

#### Simulate a single Poisson distributed trait with two variance components

$\mu = (0.2X_1 + 2.0) + X u + \epsilon$, $u \sim N(0, 0.04K)$, $\epsilon \sim N(0, 0.8I)$

$y \sim \text{Pois}(\mu)$


```julia
μ = Y ~ 0.2X1+2.0
K = cor(data')
I = eye(npeople)
Σ = [VarianceComponent(0.2, K), VarianceComponent(0.8, I)]
model = Model(μ, Σ, LogLink(), PoissonResponse())
simulate(model, data_frame)
```




<table class="data-frame"><tr><th></th><th>Y</th></tr><tr><th>1</th><td>3</td></tr><tr><th>2</th><td>3</td></tr><tr><th>3</th><td>7</td></tr><tr><th>4</th><td>0</td></tr><tr><th>5</th><td>4</td></tr><tr><th>6</th><td>1</td></tr><tr><th>7</th><td>0</td></tr><tr><th>8</th><td>21</td></tr><tr><th>9</th><td>5</td></tr><tr><th>10</th><td>7</td></tr></table>



#### A simple way to expression variance component model

Using the macro ```@vc``` instead of ```[VarianceComponent(0.2, K), VarianceComponent(0.8, I)]```


```julia
μ = Y ~ 0.2X1+2.0
K = cor(data')
I = eye(npeople)
model = Model(μ, (@vc 0.2K + 0.8I), LogLink(), PoissonResponse())
simulate(model, data_frame)
```




<table class="data-frame"><tr><th></th><th>Y</th></tr><tr><th>1</th><td>13</td></tr><tr><th>2</th><td>56</td></tr><tr><th>3</th><td>0</td></tr><tr><th>4</th><td>11</td></tr><tr><th>5</th><td>28</td></tr><tr><th>6</th><td>8</td></tr><tr><th>7</th><td>2</td></tr><tr><th>8</th><td>24</td></tr><tr><th>9</th><td>1</td></tr><tr><th>10</th><td>5</td></tr></table>



#### Simulate two traits with two variance components and cross covariance


```julia
A = [0.2 -0.1; -0.1 0.3]
B = [0.8 -0.2; -0.2 0.7]
μ = [Y1 ~ X1+0.2X2*X3+1.0, Y2 ~ X3+0.1log(HDL+LDL)+0.1]
model = Model(μ, (@vc A ⊗ K + B ⊗ I), IdentityLink(), NormalResponse(1.0))
simulate(model, data_frame)
```




<table class="data-frame"><tr><th></th><th>Y1</th><th>Y2</th></tr><tr><th>1</th><td>-0.10222997669454736</td><td>1.242014570782589</td></tr><tr><th>2</th><td>1.6985645798391524</td><td>-1.7314710486455063</td></tr><tr><th>3</th><td>0.6579712947986169</td><td>2.2033770111233073</td></tr><tr><th>4</th><td>1.0769182112301432</td><td>1.7076290809531234</td></tr><tr><th>5</th><td>2.9038228662096914</td><td>3.303032046765953</td></tr><tr><th>6</th><td>1.468375569093222</td><td>-0.31254353576852456</td></tr><tr><th>7</th><td>-1.3096109093801487</td><td>2.193982272458793</td></tr><tr><th>8</th><td>3.753520436943161</td><td>-2.1563700974996456</td></tr><tr><th>9</th><td>1.1805594385983615</td><td>0.29668922663399433</td></tr><tr><th>10</th><td>2.7257943887238394</td><td>-0.04193614515632327</td></tr></table>
