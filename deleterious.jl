using Statistics, Plots, Distributions, Random
#using Distributions, Plots, StatsPlots

"""
An individual organism
"""
struct Individual{D,L,F,N}
    deleteriousness::D
    loci::L
    fitness::F
    ndeleterious::N
end
Individual(deleteriousness::Vector, loci::Vector) = begin
    fitness = 1.0
    ndeleterious = 0
    for (i, locus) in enumerate(loci)
        if isdeleterious(locus)
            fitness *= deleteriousness[i]
            ndeleterious += 1
        end
    end
    Individual(deleteriousness, loci, fitness, ndeleterious)
end
Individual(deleteriousness::Vector, prop_deleterious::Number) = begin
    loci = [init_locus(prop_deleterious) for i in deleteriousness]
    Individual(deleteriousness, loci)
end

"""
A discrete generation of individuals.

There is no crossover between generations.
"""
struct Generation{I,O,S,F,D}
    individuals::I
    offstring::O
    size::S
    fitness::F
    ndeleterious::D
end
# Initialise from deleteriousness: individuals are sampled
Generation(deleteriousness::Vector{<:AbstractFloat}, prop_deleterious::Real, size::Int, maxsize::Int) = begin
    individuals = [Individual(deleteriousness, prop_deleterious) for i in 1:maxsize]
    offstring = [Individual(deleteriousness, prop_deleterious) for i in 1:maxsize]
    Generation(individuals, offstring, size)
end
# Initialise new generation with existing individuals
Generation(individuals::Vector{<:Individual}, offstring::Vector{<:Individual}, size::Int) = begin
    fitness = meanfitness(individuals, size)
    ndeleterious = meandeleterious(individuals, size)
    Generation(individuals, offstring, size, fitness, ndeleterious)
end

"Check if a locus is deleterious"
isdeleterious(locus) = locus[1] & locus[2]

"Choose if an allele has the deleterious gene"
init_locus(prob) = init_allele(prob), init_allele(prob)

"Choose the combination of deleterious and normal alleles at a locus"
init_allele(prob) = rand() < prob

"Randomly select a parent for an allele, from two possible"
rand_parent(a, b) = (a, b)[rand(Bool) + 1]
"Randomly slend two loci from two parents"
blend_locus(a, b) = rand_parent(a, b)[1], rand_parent(a, b)[2]
"Mate two individuals to produce another with a mixed combination of alleles"
mate!(new::Individual, a::Individual, b::Individual) = begin 
    # Reuse vectors for next generation
    deleteriousness, loci = new.deleteriousness, new.loci
    loci .= blend_locus.(a.loci, b.loci)
    Individual(deleteriousness, loci)
end

"""
Randomise the number of offspring a mother will produce, based on fitness 
and population growth rate
"""
noffspring(mother, reprorate) = rand(Poisson(reprorate * mother.fitness))

"Generate a vector of deleterious fitness values for a population"
generate_loci(nloci) = rand(0.98:1e-5:1.0, nloci)

"Run a generation"
generation!(gen::Generation, reprorate, initsize, g) = begin
    popsize = trunc(Int, initsize * reprorate^(g))
    individuals = gen.individuals
    offstring = gen.offstring
    i = 0
    tries = 0
    while i < popsize
        mother = individuals[rand(1:popsize)] 
        father = individuals[rand(1:popsize)]
        n = noffspring(mother, reprorate)
        for o in 1:n
            i += 1
            i > popsize && break
            offstring[i] = mate!(offstring[i], mother, father)
        end
        if n == 0
            tries += 1
            # Escape non-reproducing loops after 100 attemps
            if i == 0 && tries > 100
                error("Not reproducing: reduce required fitness.")
            end
        end
    end
    # Individuals vectors are recycled 
    return Generation(offstring, individuals, popsize)
end

"Calculate the mean fitness of a generation"
meanfitness(gen::Generation) = meanfitness(gen.individuals, gen.size)
meanfitness(individuals, size) = mean(individuals[i].fitness for i in 1:size)

"Calculate the mean number of deleterious loci for a generation"
meandeleterious(gen::Generation) = meandeleterious(gen.individuals, gen.size)
meandeleterious(individuals, size) = mean(individuals[i].ndeleterious for i in 1:size)

"Run a multi-generation simulation"
sim(gen::Generation, ngens::Int, reprorate::Number, nimports::Int, deleterious_loci::Vector, prop_deleterious) = begin
    meanfit = zeros(2ngens)
    meandel = zeros(2ngens)
    popsize = zeros(2ngens)
    initsize = gen.size
    initrate = 1.0
    # Run some simulations with zero growth rate
    for g = 1:ngens 
        gen = generation!(gen, initrate, initsize, g) 
        meanfit[g] = gen.fitness
        meandel[g] = gen.ndeleterious
        popsize[g] = gen.size
    end 
    # Replace some of the opulation with external individuals
    for g = gen.size+1:gen.size+nimports
        gen.individuals[g] = Individual(deleterious_loci, prop_deleterious)
    end
    gen = Generation(gen.individuals, gen.offstring, gen.size + nimports)
    # Run more simulations, now using the growth rate
    for g = 1:ngens 
        gen = generation!(gen, reprorate, initsize + nimports, g) 
        meanfit[g+ngens] = gen.fitness
        meandel[g+ngens] = gen.ndeleterious
        popsize[g+ngens] = gen.size
    end 
    return (fitness=meanfit, ndeleterious=meandel, popsize=popsize)
end

reps() = begin
    popsize = 50
    ngens = 500
    nreps = 10
    nloci = 1000
    prop_deleterious = 0.2
    nimports = 10
    reprorate = 1.01
    sims = []
    Threads.@threads for i in 1:nreps
        println("rep: ", i)
        maxsize = ceil(Int, (popsize + nimports) * reprorate^(ngens))
        deleteriousness = generate_loci(nloci)
        gen = Generation(deleteriousness, prop_deleterious, popsize, maxsize)
        push!(sims, sim(gen, ngens, reprorate, nimports, deleteriousness, prop_deleterious))
    end
    sims
end

# using Profile, ProfileView
# Profile.clear()
sims = reps()
# ProfileView.view()

ndel = [s[:ndeleterious] for s in sims]
fitness = [s[:fitness] for s in sims]
pop = sims[1].popsize

theme(:solarized)
theme(:gruvbox_light)
theme(:dark)
fitkwargs = (legend=false, xlabel="Generation", ylabel="Fitness")
delkwargs = (legend=false, xlabel="Generation", ylabel="Num deleterious loci")

plot(fitness; fitkwargs...)
plot(ndel; delkwargs...)
plot(mean(cat(fitness...; dims=2), dims=2); fitkwargs...)
plot(mean(cat(ndel...; dims=2), dims=2); delkwargs...)


