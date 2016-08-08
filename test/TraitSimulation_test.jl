module TraitSimulationTest

using TraitSimulation
using SnpArrays
using Distributions

if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
end

info("Test TraitSimulation implementation")

# load genotype data in PLINK format
genotypes = SnpData("./chr1")
num_people = genotypes.people
num_snps = genotypes.snps
snp_ids = genotypes.snpid

# choose 1000 causal SNPs at random for trait 1
srand(0)
num_causals = 1000
causals = snp_ids[sort(randperm(num_snps)[1:num_causals])]
hsq1 = 0.2
normal_dist = Normal(0.0, sqrt(hsq1 / num_causals))
betas = rand(normal_dist, num_causals)
causal_betas = Dict{AbstractString, Float64}()
for i in 1:num_causals
  causal_betas[causals[i]] = betas[i]
end

# choose 1000 causal SNPs at random for trait 2
num_causals = 2000
causals = snp_ids[sort(randperm(num_snps)[1:num_causals])]
hsq2 = 0.5
normal_dist = Normal(0.0, sqrt(hsq2 / num_causals))
gammas = rand(normal_dist, num_causals)
causal_gammas = Dict{AbstractString, Float64}()
for i in 1:num_causals
  causal_gammas[causals[i]] = gammas[i]
end

# simulate a single quantitative trait
simulate(genotypes, :q, causal_betas, hsq1, rep=5, missing_rate=0.01)

# simulate a single case-control trait
simulate(genotypes, :b, causal_betas, hsq1, rep=5, prevalence=0.25,
  ncc=(50, 100))

# simulate two quantitative traits
simulate(genotypes, (:q,:q), (causal_betas, causal_gammas),
  (hsq1,hsq2), trait_cor=0.8, rep=5, missing_rate=(0.01,0.02))

# simulate a quantitative trait and a binary trait
simulate(genotypes, (:q,:b), (causal_betas, causal_gammas),
  (hsq1,hsq2), trait_cor=0.2, rep=5, prevalence=0.25,
  missing_rate=0.01, ncc=(50, 100))

# simulate a case-control trait and a quantitative trait
simulate(genotypes, (:b,:q), (causal_betas, causal_gammas),
    (hsq1,hsq2), trait_cor=0.2, rep=5, prevalence=0.25,
    missing_rate=0.01, ncc=(50, 100))

# simulate two case-control traits
simulate(genotypes, (:b,:b), (causal_betas, causal_gammas),
    (hsq1,hsq2), trait_cor=0.2, rep=5, prevalence=(0.25, 0.3),
    ncc=((50, 100),(60, 110)))

end # end module
