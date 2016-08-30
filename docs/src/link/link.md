# Link Functions

This section lists the link function types supported by the TraitSimulation
module.

## LinkFunction

```LinkFunction``` is an abstract type. It's the super type of each
specific link function type listed below. All sub types of ```LinkFunction```
have constructors of the form ```[name of the link]()```, and
contain ```link_inv::Function``` as their member variables. For example
the constructor for identity link function type is ```IdentityLink()```,
and its ```link_inv``` member variable is the function ```f(x) = x```.

## LinkFunctionType

```LinkFunctionType``` is a type alias for a single link function or a
vector of link functions.

## CauchitLink

Implements the Cauchit link $g(\mu) = \text{arctan}(\mu)/\pi+1/2$,
which has inverse $g^{-1}(\eta) = \text{tan}[(\eta-1/2)\pi]$.

## CloglogLink

Implements the complementary log log link, $g(\mu) = \log(-\log(1-\mu))$
which has invserse $g^{-1}(\eta) = 1 - \exp(-\exp(\eta))$.

## IdentityLink

Implements the identity link $g(\mu) = \mu$, which has inverse
$g^{-1}(\eta) = \eta$

## InverseLink

Implements the inverse link $g(\mu) = 1/\mu$, which has inverse
$g^{-1}(\eta) = 1/\eta$

## LogitLink

Implements the logit link $g(\mu) = \log(\mu / (1-\mu))$, which has inverse
$g^{-1}(\eta) = 1 / (1+\exp(-\eta))$.

## ProbitLink

Implements the probit link $g(\mu) = \Phi^{-1}(\mu)$,
which has inverse $g^{-1}(\eta) = \Phi(\eta)$

## SqrtLink

Implements the square root link $g(\mu) = \sqrt{\mu}$, which has inverse
$g^{-1}(\eta) = \eta^2$.

## LogLink

Implements the log link $g(\mu) = \log(\mu)$, which has inverse
$g^{-1}(\eta) = \exp(\eta)$.


