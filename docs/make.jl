using Documenter, TraitSimulation

makedocs()

deploydocs(
  deps   = Deps.pip("pygments", "mkdocs", "mkdocs-material", "python-markdown-math"),
  repo   = "https://github.com/huwenboshi/TraitSimulation.jl.git",
  julia  = "0.5",
  osname = "linux"
  )
