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


####Create a random data set for testing

```julia
using DataFrames, TraitSimulation
# create a dataset for testing
srand(1)
data = rand(5,6)
df = convert(DataFrame, data)
cov_mat = cor(data')
names!(df, [:A, :B, :C, :D, :E, :F])
```




<table class="data-frame"><tr><th></th><th>A</th><th>B</th><th>C</th><th>D</th><th>E</th><th>F</th></tr><tr><th>1</th><td>0.23603334566204692</td><td>0.21096820215853596</td><td>0.5557510873245723</td><td>0.20947237319807077</td><td>0.07695088688120899</td><td>0.6448833539420931</td></tr><tr><th>2</th><td>0.34651701419196046</td><td>0.951916339835734</td><td>0.43710797460962514</td><td>0.25137920979222494</td><td>0.6403962459899388</td><td>0.07782644396003469</td></tr><tr><th>3</th><td>0.3127069683360675</td><td>0.9999046588986136</td><td>0.42471785049513144</td><td>0.02037486871266725</td><td>0.8735441302706854</td><td>0.8481854810000327</td></tr><tr><th>4</th><td>0.00790928339056074</td><td>0.25166218303197185</td><td>0.773223048457377</td><td>0.2877015122756894</td><td>0.27858242002877853</td><td>0.0856351682044918</td></tr><tr><th>5</th><td>0.4886128300795012</td><td>0.9866663668987996</td><td>0.2811902322857298</td><td>0.859512136087661</td><td>0.7513126327861701</td><td>0.5532055454580578</td></tr></table>



####Simulate a single trait with identity link and normal response
```julia
# simulate a single trait using GLM
sim_model = Model(T ~ A + 2B*C, IdentityLink(), NormalResponse(1.0))
y = simulate(sim_model, df)
```




<table class="data-frame"><tr><th></th><th>T</th></tr><tr><th>1</th><td>0.8296713687454297</td></tr><tr><th>2</th><td>1.4987221765068286</td></tr><tr><th>3</th><td>1.4212780112467271</td></tr><tr><th>4</th><td>0.856787269123028</td></tr><tr><th>5</th><td>1.2534810799152252</td></tr></table>



####Simulate two traits with the same link and response distribution
```julia
# simulate two traits with the same link and response using GLM
formulae = [T1 ~ A+2B*C, T2 ~ A + 2log(B+C)+2.0]
sim_model = Model(formulae, IdentityLink(), NormalResponse(1.0))
y = simulate(sim_model, df)
```




<table class="data-frame"><tr><th></th><th>T1</th><th>T2</th></tr><tr><th>1</th><td>1.9937003523203445</td><td>-0.005103919607121421</td></tr><tr><th>2</th><td>2.4384920644109673</td><td>4.5036580275445015</td></tr><tr><th>3</th><td>0.9967785820767805</td><td>2.8916050756826595</td></tr><tr><th>4</th><td>0.5532857789690359</td><td>1.6068865859160373</td></tr><tr><th>5</th><td>0.2291794254963837</td><td>1.4165366031669122</td></tr></table>



####Simulate three traits with different links and response distributions
```julia
# simulate three traits with different link and response using GLM
formulae = [T1 ~ A+2B*C, T2 ~ A+2log(B+C)+2.0, T3 ~ A+B+C+1.0]
links = [IdentityLink(), LogitLink(), LogLink()]
resp_dists = [NormalResponse(1.0), BinomialResponse(100), PoissonResponse()]
sim_model = Model(formulae, links, resp_dists)
y = simulate(sim_model, df)
```




<table class="data-frame"><tr><th></th><th>T1</th><th>T2</th><th>T3</th></tr><tr><th>1</th><td>-1.4962208991713917</td><td>88</td><td>11</td></tr><tr><th>2</th><td>2.603781693021933</td><td>96</td><td>16</td></tr><tr><th>3</th><td>0.7634421808186758</td><td>98</td><td>14</td></tr><tr><th>4</th><td>1.2624503185897922</td><td>87</td><td>10</td></tr><tr><th>5</th><td>-0.48139210111585573</td><td>96</td><td>23</td></tr></table>



####Simulate a single trait under GLMM with two variance components
```julia
# simulate a single trait with two variance components using GLMM
vc = [VarianceComponent(0.2, cov_mat), VarianceComponent(0.8, eye(5))]
sim_model = Model(T ~ A+2B*C, LogitLink(), BinomialResponse(100), vc)
y = simulate(sim_model, df)
```




<table class="data-frame"><tr><th></th><th>T</th></tr><tr><th>1</th><td>68</td></tr><tr><th>2</th><td>62</td></tr><tr><th>3</th><td>76</td></tr><tr><th>4</th><td>43</td></tr><tr><th>5</th><td>72</td></tr></table>



####Simulate two traits each with two variance components
```julia
# simulate two traits with two variance components using GLMM
# the two traits can have different response distribution
formulae = [T1 ~ A+2B*C, T2 ~ C+log(C)+3.0]
links = [IdentityLink(), LogitLink()]
resp_dists = [NormalResponse(1.0), PoissonResponse()]
vc = [VarianceComponent([0.2, 0.3], cov_mat),
      VarianceComponent([0.8, 0.7], eye(5))]
sim_model = Model(formulae, links, resp_dists, vc)
y = simulate(sim_model, df)
```




<table class="data-frame"><tr><th></th><th>T1</th><th>T2</th></tr><tr><th>1</th><td>2.300793397259024</td><td>1</td></tr><tr><th>2</th><td>-0.9724957256066538</td><td>0</td></tr><tr><th>3</th><td>0.16186818794850574</td><td>0</td></tr><tr><th>4</th><td>-1.981011982636921</td><td>1</td></tr><tr><th>5</th><td>-1.5567734788160243</td><td>0</td></tr></table>



####Simulate two traits with cross covariances
```julia
# simulate two traits with two variance components with cross covariances
# using GLMM. the two traits have the same response distribution
formulae = [T1 ~ A+2B*C, T2 ~ C+log(C)+3.0]
vc = [VarianceComponent([0.2 -0.1; -0.1 0.3], cov_mat),
      VarianceComponent([0.8 -0.2; -0.2 0.7], eye(5))]
sim_model = Model(formulae, IdentityLink(), NormalResponse(1.0), vc)
y = simulate(sim_model, df)
```




<table class="data-frame"><tr><th></th><th>T1</th><th>T2</th></tr><tr><th>1</th><td>1.2375761672364602</td><td>5.598899068351228</td></tr><tr><th>2</th><td>1.7749641269608285</td><td>0.4235391538698168</td></tr><tr><th>3</th><td>2.0820850495371017</td><td>0.20826081606112956</td></tr><tr><th>4</th><td>-0.41284142811269164</td><td>3.5844275050705</td></tr><tr><th>5</th><td>0.3612051541311363</td><td>0.5470808314522484</td></tr></table>
