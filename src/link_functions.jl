"""
A list of types to store inverse link functions
"""


# super type for all link function types
abstract LinkFunction


# Cauchit link function
type CauchitLink <: LinkFunction
  link::Function
  link_inv::Function
end

CauchitLink() = CauchitLink(x::Float64 -> atan(x)/pi + 0.5,
                            x::Float64 -> tan(pi*(x-0.5)))


# complementary log log link function
type CloglogLink <: LinkFunction
  link::Function
  link_inv::Function
end

CloglogLink() = CloglogLink(x::Float64 -> log(-log(1-x)),
                            x::Float64 -> 1.0-exp(-exp(x)))


# identity link function
type IdentityLink <: LinkFunction
  link::Function
  link_inv::Function
end

IdentityLink() = IdentityLink(x::Float64 -> x,
                              x::Float64 -> x)


# inverse link function
type InverseLink <: LinkFunction
  link::Function
  link_inv::Function
end

InverseLink() = InverseLink(x::Float64 -> 1.0/x,
                            x::Float64 -> 1.0/x)


# logit link function
type LogitLink <: LinkFunction
  link::Function
  link_inv::Function
end

LogitLink() = LogitLink(x::Float64 -> log(x/(1.0-x)),
                        x::Float64 -> 1.0/(1.0+exp(-x)))


# probit link function
const normal_dist = Normal(0.0, 1.0)
type ProbitLink <: LinkFunction
  link::Function
  link_inv::Function
end

ProbitLink() = ProbitLink(x::Float64 -> quantile(normal_dist, x),
                          x::Float64 -> cdf(normal_dist, x))


# square root link function
type SqrtLink <: LinkFunction
  link::Function
  link_inv::Function
end

SqrtLink() = SqrtLink(x::Float64 -> √x, 
                      x::Float64 -> x*x)


# log link function
type LogLink <: LinkFunction
  link::Function
  link_inv::Function
end

LogLink() = LogLink(x::Float64 -> log(x),
                    x::Float64 -> exp(x))
