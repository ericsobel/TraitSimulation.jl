

```julia
include("../src/TraitSimulation.jl")
using DataFrames, TraitSimulation
```


```julia
df = convert(DataFrame, rand(10,6))
names!(df, [:A, :B, :C, :D, :E, :F])
```




<table class="data-frame"><tr><th></th><th>A</th><th>B</th><th>C</th><th>D</th><th>E</th><th>F</th></tr><tr><th>1</th><td>0.4933257820869694</td><td>0.8739836725086518</td><td>0.5781556567086412</td><td>0.13411863983952643</td><td>0.3037671612591084</td><td>0.6779496718496925</td></tr><tr><th>2</th><td>0.40230119346414783</td><td>0.5666733787500011</td><td>0.6871510261988054</td><td>0.9215853800851344</td><td>0.9096225210283237</td><td>0.9679497165073165</td></tr><tr><th>3</th><td>0.8899053655154903</td><td>0.09647328170995118</td><td>0.623114899405943</td><td>0.08660519400659727</td><td>0.15880266674611265</td><td>0.05448046581645949</td></tr><tr><th>4</th><td>0.43479090551664124</td><td>0.6627991742690582</td><td>0.5784006898435621</td><td>0.24297453610281683</td><td>0.21490937224447437</td><td>0.2569346597765203</td></tr><tr><th>5</th><td>0.8617670737428211</td><td>0.1904768111531845</td><td>0.4825223264952174</td><td>0.9657665868006096</td><td>0.7181465303396555</td><td>0.9753196752861237</td></tr><tr><th>6</th><td>0.9175237947348043</td><td>0.09087254000144385</td><td>0.6015962171510085</td><td>0.7320931219503042</td><td>0.882685513480935</td><td>0.9513196544496756</td></tr><tr><th>7</th><td>0.10520195436570279</td><td>0.9798124633577467</td><td>0.2722336965394485</td><td>0.3053081097480914</td><td>0.39163865020696487</td><td>0.19695641970029532</td></tr><tr><th>8</th><td>0.562779752014261</td><td>0.3511162863985213</td><td>0.11062369559289609</td><td>0.9606634940577203</td><td>0.8829051340487035</td><td>0.25146286321576516</td></tr><tr><th>9</th><td>0.6238790308934099</td><td>0.2651779533764913</td><td>0.567448521248149</td><td>0.7515804495315923</td><td>0.2018242811008908</td><td>0.6418052416557587</td></tr><tr><th>10</th><td>0.38528844811331453</td><td>0.5366003551865897</td><td>0.19310447368685524</td><td>0.8299717338469788</td><td>0.5408173192325361</td><td>0.6840383010171518</td></tr></table>




```julia
formula = T ~ A+2B*C+log(3D*(E+0.8F))+2.0
```




    Formula: T ~ A + (2B) * C + log((3D) * (E + 0.8F)) + 2.0




```julia
# simulate a normal response with σ=1.0
sim_model = Model(formula, IdentityLink(), NormalResponse(1.0))
simulate(sim_model, df)
```




<table class="data-frame"><tr><th></th><th>T</th></tr><tr><th>1</th><td>3.504505164086454</td></tr><tr><th>2</th><td>4.310869399429139</td></tr><tr><th>3</th><td>1.0876713317408702</td></tr><tr><th>4</th><td>2.3914199497591464</td></tr><tr><th>5</th><td>4.4325351124422445</td></tr><tr><th>6</th><td>2.8993251075190107</td></tr><tr><th>7</th><td>2.673017795010866</td></tr><tr><th>8</th><td>4.021456711699673</td></tr><tr><th>9</th><td>3.4193052423320878</td></tr><tr><th>10</th><td>5.17637945995563</td></tr></table>




```julia
# simulate a binomial response with n=100
sim_model = Model(formula, LogitLink(), BinomialResponse(100))
simulate(sim_model, df)
```




<table class="data-frame"><tr><th></th><th>T</th></tr><tr><th>1</th><td>93</td></tr><tr><th>2</th><td>99</td></tr><tr><th>3</th><td>49</td></tr><tr><th>4</th><td>91</td></tr><tr><th>5</th><td>100</td></tr><tr><th>6</th><td>99</td></tr><tr><th>7</th><td>88</td></tr><tr><th>8</th><td>97</td></tr><tr><th>9</th><td>98</td></tr><tr><th>10</th><td>96</td></tr></table>




```julia
# simulate a Poisson response
sim_model = Model(formula, LogLink(), PoissonResponse())
simulate(sim_model, df)
```




<table class="data-frame"><tr><th></th><th>T</th></tr><tr><th>1</th><td>12</td></tr><tr><th>2</th><td>99</td></tr><tr><th>3</th><td>1</td></tr><tr><th>4</th><td>17</td></tr><tr><th>5</th><td>99</td></tr><tr><th>6</th><td>77</td></tr><tr><th>7</th><td>10</td></tr><tr><th>8</th><td>37</td></tr><tr><th>9</th><td>32</td></tr><tr><th>10</th><td>33</td></tr></table>




```julia
# simulate a Bernoulli response
sim_model = Model(formula, LogitLink(), BernoulliResponse())
simulate(sim_model, df)
```




<table class="data-frame"><tr><th></th><th>T</th></tr><tr><th>1</th><td>1</td></tr><tr><th>2</th><td>1</td></tr><tr><th>3</th><td>0</td></tr><tr><th>4</th><td>1</td></tr><tr><th>5</th><td>1</td></tr><tr><th>6</th><td>1</td></tr><tr><th>7</th><td>1</td></tr><tr><th>8</th><td>1</td></tr><tr><th>9</th><td>1</td></tr><tr><th>10</th><td>1</td></tr></table>




```julia
# simulate an Exponential response
sim_model = Model(formula, InverseLink(), ExponentialResponse())
simulate(sim_model, df)
```




<table class="data-frame"><tr><th></th><th>T</th></tr><tr><th>1</th><td>0.02375534822371872</td></tr><tr><th>2</th><td>0.14055072607201588</td></tr><tr><th>3</th><td>8.876709263297036</td></tr><tr><th>4</th><td>0.3468466653555284</td></tr><tr><th>5</th><td>0.35110397796023435</td></tr><tr><th>6</th><td>0.30856165093118143</td></tr><tr><th>7</th><td>0.08373521564397778</td></tr><tr><th>8</th><td>0.2534106842926165</td></tr><tr><th>9</th><td>0.4525316281929459</td></tr><tr><th>10</th><td>0.026888275196969787</td></tr></table>




```julia
# simulate a gamma response with shape parameter 2.0
sim_model = Model(formula, InverseLink(), GammaResponse(2.0))
simulate(sim_model, df)
```




<table class="data-frame"><tr><th></th><th>T</th></tr><tr><th>1</th><td>0.34434393232264326</td></tr><tr><th>2</th><td>0.2890975230054341</td></tr><tr><th>3</th><td>14.830712722934773</td></tr><tr><th>4</th><td>0.9323246533111741</td></tr><tr><th>5</th><td>0.6516330211627606</td></tr><tr><th>6</th><td>0.460886898208853</td></tr><tr><th>7</th><td>1.2446424985637627</td></tr><tr><th>8</th><td>1.5788506659265176</td></tr><tr><th>9</th><td>0.1967905622466749</td></tr><tr><th>10</th><td>0.13124277622658273</td></tr></table>




```julia

```
