module GWASSimulationTestQuantitativeSingle

using GWASSimulation
using SnpArrays
using Distributions

if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
end

info("Test GWASSimulation implementation - a pair of quantitative traits")

# load genotype data in PLINK format
genotypes = SnpData("./chr1")
num_people = genotypes.people
num_snps = genotypes.snps
snp_ids = genotypes.snpid

# choose 1000 causal SNPs at random for trait 1
srand(0)
num_causals = 1000
causals = snp_ids[sort(randperm(num_snps)[1:num_causals])]

# sample effect size beta from normal distribution for trait 1
hsq1 = 0.2
normal_dist = Normal(0.0, sqrt(hsq1 / num_causals))
betas = rand(normal_dist, num_causals)
causal_betas = Dict{AbstractString, Float64}()
for i in 1:num_causals
  causal_betas[causals[i]] = betas[i]
end

# choose 1000 causal SNPs at random for trait 2
srand(0)
num_causals = 2000
causals = snp_ids[sort(randperm(num_snps)[1:num_causals])]

# sample effect size beta from normal distribution for trait 2
hsq2 = 0.5
normal_dist = Normal(0.0, sqrt(hsq2 / num_causals))
gammas = rand(normal_dist, num_causals)
causal_gammas = Dict{AbstractString, Float64}()
for i in 1:num_causals
  causal_gammas[causals[i]] = gammas[i]
end


# number of replications to simulate
num_reps = 5
missing_rate1 = 0.01
missing_rate2 = 0.02

# simulate a quantitative trait
simulate(genotypes, (:q,:q), (causal_betas, causal_gammas),
  (hsq1,hsq2), trait_cor=0.2, rep=num_reps,
  missing_rate=(missing_rate1,missing_rate2))

end # end module
