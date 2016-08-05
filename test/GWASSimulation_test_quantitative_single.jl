module GWASSimulationTestQuantitativeSingle

using GWASSimulation
using SnpArrays
using Distributions

if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
end

info("Test GWASSimulation implementation")

# load genotype data in PLINK format
genotypes = SnpData("./chr1")
num_people = genotypes.people
num_snps = genotypes.snps
snp_ids = genotypes.snpid

# choose 1000 causal SNPs at random
srand(0)
num_causals = 1000
causals = snp_ids[sort(randperm(num_snps)[1:num_causals])]

# sample effect size beta from normal distribution
hsq = 0.2
normal_dist = Normal(0.0, sqrt(hsq / num_causals))
betas = rand(normal_dist, num_causals)
causal_betas = Dict{AbstractString, Float64}()
for i in 1:num_causals
  causal_betas[causals[i]] = betas[i]
end

# number of replications to simulate
num_reps = 5
missing_rate = 0.01

# simulate a quantitative trait
simulate(genotypes, :q, causal_betas, hsq, rep=num_reps,
  missing_rate=missing_rate)

end # end module
