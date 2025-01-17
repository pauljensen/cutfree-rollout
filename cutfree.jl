

import Base: *, sort
using Random, Statistics, DataFrames, CSV
using BioSequences

# Overload the concatenation operator to combine a sequence and 
# a single base, i.e. dna"AGCGTGC" * DNA_T
function *(seq::LongDNASeq, base::DNA)
    seq * LongDNASeq([base])
end

# Overload sort to order DNA codes by degeneracy:
#   A,G,C,T,M,R,W,S,Y,K,V,H,D,B,N
function sort(seq::LongSequence{DNAAlphabet{4}}; rev::Bool=false)
    sort(convert(Vector{DNA}, seq), by=degeneracy, rev=rev)
end

"""
    isvalid(seq, sites)

Check if any of the sites appear in the sequence.
"""
function isvalid(seq, sites)
    for site in sites
        if occursin(site, seq)
            return false
        end
    end
    return true
end

"""
    find_valid_randomer(prefix, bases, sites)

Randomly select bases to add the prefix, ensuring no site
appears in the randomer.
"""
function find_valid_randomer(prefix, bases, sites)
    n = length(bases)
    pre_n = length(prefix)
    randomer = prefix * (dna"-" ^ n)
    for i = 1:n
        candidates = shuffle(bases[i])
        for cand in candidates
            randomer[pre_n+i] = cand
            if isvalid(randomer[1:pre_n+i], sites)
                break
            else
                randomer[pre_n+i] = DNA_Gap
            end
        end
        if randomer[pre_n+i] == DNA_Gap
            break
        end
    end
    return randomer
end

"""
    degeneracy(base::DNA[, uselog2=true])

Calculate the degeneracy of a single base. The degeneracy 
is the number of bases in the code, so
    degeneracy(DNA_N) == 4
and
    degeneracy(DNA_A) == 1

If `uselog2`, the degeneracy is log2 transformed. By convention,
the degneracy of `DNA_Gap` is zero in both the log2-transformed
and the untransformed case.
"""
function degeneracy(base::DNA; uselog2=true)
    deg = 0.0
    if base == DNA_N
        deg = 4.0
    elseif base == DNA_A || base == DNA_C || base == DNA_G || base == DNA_T
        deg = 1.0
    elseif base == DNA_B || base == DNA_D || base == DNA_H || base == DNA_V
        deg = 3.0
    elseif (base == DNA_K || base == DNA_Y || base == DNA_S || base == DNA_W 
            || base == DNA_R || base == DNA_M)
        deg = 2.0
    else # DNA_Gap
        deg = 0.0
    end
    
    if uselog2
        if base == DNA_Gap
            deg = 1.0
        end
        return log2(deg)
    else
        return deg
    end
end

"""
    degeneracy(seq::LongSequence{DNAAlphabet{4}}; uselog2=true)

Calculate the degeneracy of a DNA sequence.
"""
function degeneracy(seq::LongSequence{DNAAlphabet{4}}; uselog2=true)
    if uselog2
        return sum(degeneracy.(seq, uselog2=true))
    else
        return prod(degeneracy.(seq))
    end
end

# Run Monte Carlo simulations using a random base policy over a fixed horizon.
function simulate_random(prefix, bases, sites; nsims=1000, horizon=length(bases))
    degs = zeros(nsims)
    true_horizon = min(horizon, length(bases)) # don't go beyond the randomer length
    for i = 1:nsims
        degs[i] = degeneracy(find_valid_randomer(prefix, bases[1:true_horizon], sites))
    end
    return mean(degs)
end

# Run a 1-step greedy lookahead policy.
function simulate_greedy(prefix, bases, sites; horizon=length(bases))
    true_horizon = min(horizon, length(bases)) # don't go beyond the randomer length
    for h = 1:true_horizon
        # Iterate through the bases by decending degeneracy, so the greedy approach
        # is to stop as soon as we've found one.
        found = false # define out here for scoping rules
        for base in sort(bases[h], rev=true)
            found = false
            if isvalid(prefix * base, sites)
                prefix *= base
                found = true
                break
            end
        end
        if ~found
            # None of the bases were valid, so stop and leave the DNA_Gap
            break
        end
    end
    return degeneracy(prefix)
end


function cutfree_rollout(bases, sites; simulate=simulate_random, kwargs...)
    n = length(bases)
    randomer = dna"-" ^ n

    # If the restriction site is non-palindromic, we also need to block
    # the reverse complement.
    for i = 1:length(sites)
        if ~ispalindromic(sites[i])
            push!(sites, reverse_complement(sites[i]))
        end
    end

    # Find the longest restriction site and use this as the horizon.
    max_len = map(length, sites) |> maximum
    horizon = max_len

    for i = 1:n
        best_base = DNA_Gap
        best_deg = -1.0
        start = max(i - max_len, 1) # only need to look back max_len when checking the oligos
        for base in bases[i]
            if ~isvalid(randomer[start:i-1] * base, sites)
                # adding this bases creates a restriction site; skip it
                continue
            end
            deg = simulate(randomer[start:i-1] * base, bases[i+1:end], sites; horizon=horizon, kwargs...)

            # Update use this base if:
            #   1. the mean randomer degeneracy is higher than the previous best, OR
            #   2. the randomer degeneracy is tied with the previous best, but the new base has
            #      higher individual degeneracy
            if deg > best_deg || (deg == best_deg && degeneracy(base) > degeneracy(best_base))
                best_deg = deg
                best_base = base
            end
        end

        if best_base == DNA_Gap
            # The sequence has terminated; there are no bases that can be added without making
            # a restriction site.
            break
        else
            randomer[i] = best_base
        end
    end

    return randomer
end


all_bases = dna"AGCTMRWSYKVHDBN"
all_ns = [all_bases for _ in 1:10] # make a length 10 randomer
sites = [dna"GGTCTC"] # BsaI

randomer = cutfree_rollout(all_ns, sites, simulate=simulate_random, nsims=1000)
randomer = cutfree_rollout(all_ns, sites, simulate=simulate_greedy)

deg = degeneracy(randomer)
natdeg = deg / log2(exp(1))
println("Final randomer: $randomer")
println("Final degeneracy: $deg ($natdeg)")

# Run the benchmark CutFree set
function runtimes(data)
    data[:random_objval] = 0.0
    data[:random_runtime] = 0.0
    data[:random_code] = ""

    data[:greedy_objval] = 0.0
    data[:greedy_runtime] = 0.0
    data[:greedy_code] = ""

    data.objval = data.objval ./ log(2)

    bases = [all_bases for _ in 1:20]
    sites = map(x -> split(x, ',') .|> LongDNASeq, data.sites)
    for i = 1:nrow(data)
        stats = @timed cutfree_rollout(bases, sites[i], nsims=1000)
        data[i,:random_objval] = degeneracy(stats.value)
        data[i,:random_runtime] = stats.time
        data[i,:random_code] = convert(String, stats.value)

        stats = @timed cutfree_rollout(bases, sites[i], simulate=simulate_greedy)
        data[i,:greedy_objval] = degeneracy(stats.value)
        data[i,:greedy_runtime] = stats.time
        data[i,:greedy_code] = convert(String, stats.value)
        println(i)
    end
    return data
end

#data = runtimes(CSV.read("./results.csv")[1:10,:])





