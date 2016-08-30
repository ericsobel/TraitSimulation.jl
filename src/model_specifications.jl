"""
The Model abstract type, super type for each specific simulation model
"""
abstract SimulationModel

# Define some aliases for clarity
typealias FormulaType
  Union{Formula, Vector{Formula}}

typealias VarianceComponentType
  Union{Float64, Vector{Float64}, Matrix{Float64}}

typealias TraitType
  Union{Symbol, Vector{Symbol}}

"""
A type to store fixed effect simulation model
"""
type FixedEffectModel <: SimulationModel

  # specify the formula of the simulation
  formula::FormulaType
  
  # specify the link function
  link::LinkFunctionType

  # specify the distribution of the response
  resp_dist::ResponseDistributionType

end


# return the number of traits being simulated
Base.size(model::FixedEffectModel) = 
typeof(model.formula) == Formula ? 1 : size(model.formula, 1)


"""
A type to store variance component and its covariance matrix
"""
# TODO: Implement a "toString" function for print
type VarianceComponent

  # Stores a single variance component, a vector of variance component
  #for a number of traits, or a cross covariance matrix
  var_comp::VarianceComponentType

  # stores the covariance matrix
  cov_mat::Matrix{Float64}

end


"""
A type to store random effect simulation model
"""
type RandomEffectModel <: SimulationModel

  # specify trait names
  traits::TraitType

  # specify the variance components for GLMM
  vc::Vector{VarianceComponent}
  
  # specify the link function
  link::LinkFunctionType

  # specify the distribution of the response
  resp_dist::ResponseDistributionType

end


# return the number of traits being simulated
Base.size(model::RandomEffectModel) =
typeof(model.vc[1].var_comp) == Float64 ? 1 : size(model.vc.var_comp, 1)


"""
A type to store fixed effect simulation model
"""
type MixedEffectModel <: SimulationModel

  # specify the formula of the simulation
  formula::FormulaType

  # specify the variance components for GLMM
  vc::Vector{VarianceComponent}
  
  # specify the link function
  link::LinkFunctionType

  # specify the distribution of the response
  resp_dist::ResponseDistributionType

end


# return the number of traits being simulated
Base.size(model::MixedEffectModel) = 
typeof(model.formula) == Formula ? 1 : size(model.formula, 1)


"""
Parse an expression specifying the covariance matrix of the random effects
Returns code for constructing the array of VarianceComponent
"""
macro vc(expr::Expr)
  
  # TODO: add code to check validity of expr before parseing it
  # throw an exception if the expr doesn't follow the syntax

  expr_str = string(expr)
  
  ret_expr = quote
    [VarianceComponent(eval(parse($expr_str).args[i].args[2]),
                       eval(parse($expr_str).args[i].args[3]))
     for i=2:size(parse($expr_str).args,1)]
  end

  return esc(ret_expr)

end
