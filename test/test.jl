import JuliaParser.Parser

"""
Stores model of a formula
"""
type ModelFormula
  lhs::Symbol
  rhs::Expr
end

"""
Construct a model formula following DataFrames convention, except that
coefficients are allowed in front of every term.
"""
macro ~(lhs, rhs)
  ex = Expr(:call, :ModelFormula, Base.Meta.quot(lhs), Base.Meta.quot(rhs))
  return ex
end

exp = y ~ x+1.5log(y+z)
ast = Parser.parse(exp.rhs)

println("hell world.")
