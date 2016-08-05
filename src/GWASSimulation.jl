module GWASSimulation

export Simulation, simulate

using SnpArrays
using Distributions

"""
"""
type Simulation
end

"""
"""
function simulate(
    genotypes::SnpData,
    trait_type::Union{Symbol, NTuple{2, Symbol}},
    causal_snps::Union{Dict{AbstractString, Float64},
      Dict{AbstractString, NTuple{2,Float64}}},
    hsq::Union{Float64, NTuple{2, Float64}};
    rep::Int64 = 1,
    missing_rate::Float64 = 0.0,
    prevalence::Union{Float64, NTuple{2, Float64}} = 0.0,
    cor::Float64 = 0.0,
    ncc::Union{NTuple{2, Int64}, NTuple{4, Int64}} = (0,0)
  )

  # simulate a single quantitative trait
  if (typeof(trait_type) == Symbol && trait_type == :q)
    simulate_liability(genotypes, causal_snps, hsq, rep, missing_rate)
    return
  end

  # simulate a single binary (case-control) trait
  if (typeof(trait_type) == Symbol && trait_type == :b)
    simulate_case_control(genotypes, causal_snps, hsq, rep,
      prevalence, ncc[1], ncc[2])
    return
  end

  # simulate a pair of correlated quantitative traits

  # simulate a quantitative trait and a binary (case-control) trait

  # simulate a pair of binary (case-control) traits

end # end function simulate

"""
"""
function simulate_liability(
    genotypes::SnpData,
    causal_snps::Dict{AbstractString, Float64},
    hsq::Union{Float64, NTuple{2, Float64}},
    rep::Int64,
    missing_rate::Float64
  )

  # get dimensions
  num_snps = genotypes.snps
  num_people = genotypes.people

  # get centered and scaled genotype matrix X
  X = convert(Array{Float64, 2}, genotypes.snpmatrix,
        impute=true, center=true, scale=true)

  # construct β
  β = zeros(Float64, num_snps)
  for i in 1:num_snps
    if (haskey(causal_snps, genotypes.snpid[i]))
      β[i] = causal_snps[genotypes.snpid[i]]
    end
  end

  # genetic component of the trait
  Xβ = X * β

  # variance of the genetic component of the trait
  σ2g = var(Xβ)

  # compute σ2e so that σ2g / (σ2g + σ2e) = hsq
  σ2e = σ2g * (1.0 - hsq) / hsq

  # modeling environmental component of the trait
  norm_dist = Normal(0.0, sqrt(σ2e))

  # modeling missing rate
  bern_dist = Bernoulli(missing_rate)

  # simulate replicates
  Y = zeros(Float64, (num_people, rep))
  for i in 1:rep

    # draw environmental effect and generate phenotype
    ϵ = rand(norm_dist, num_people)
    Y[:,i] = Xβ + ϵ

    # mark traits as missing based on missing_rate
    for j in 1:num_people
      if (rand(bern_dist) == 1)
        Y[j,i] = NaN
      end
    end

  end

  return Y

end # end function sim_liability

"""
"""
function simulate_case_control(
    genotypes::SnpData,
    causal_snps::Dict{AbstractString, Float64},
    hsq::Union{Float64, NTuple{2, Float64}},
    rep::Int64,
    prevalence::Float64,
    num_case::Int64,
    num_control::Int64
  )

  # get dimensions
  num_people = genotypes.people

  # simulate the liability scores
  Y_liab = simulate_liability(genotypes, causal_snps, hsq, rep, 0.0)

  # mark cases and controls
  case_status = falses((num_people, rep))
  for i in 1:rep

    # get liability threshold
    k::Int64 = ceil(prevalence * num_people)
    sorted = sort(Y_liab[:,i], rev=true)
    threshold = sorted[k]

    # mark case and control status based on threshold
    for j in 1:num_people
      if (Y_liab[j,i] > threshold)
        case_status[j,i] = true
      end
    end

  end

  # simulate ascertainment bias based on number of cases and controls
  Y = zeros(Float64, (num_people, rep))
  for i in 1:rep

    # get case and control indices
    case_indices = find(case_status[:,i] .== false)
    ctrl_indices = find(case_status[:,i] .== true)

    

  end

end



end # end module
