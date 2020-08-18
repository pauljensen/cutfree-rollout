

import Base: *
using Random, Statistics, DataFrames, CSV
using BioSequences

function is_valid(seq, sites)
    for site in sites
        if occursin(site, seq)
            return false
        end
    end
    return true
end

function find_valid_randomer(prefix, bases, sites)
    n = length(bases)
    pre_n = length(prefix)
    randomer = prefix * (dna"-" ^ n)
    for i = 1:n
        candidates = shuffle(bases[i])
        for cand in candidates
            randomer[pre_n+i] = cand
            if is_valid(randomer[1:pre_n+i], sites)
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

function degeneracy(seq)
    deg = 0.0
    for i = 1:length(seq)
        if seq[i] == DNA_N
            deg += 2.0
        elseif seq[i] == DNA_A || seq[i] == DNA_C || seq[i] == DNA_G || seq[i] == DNA_T
            deg += 0
        elseif seq[i] == DNA_B || seq[i] == DNA_D || seq[i] == DNA_H || seq[i] == DNA_V
            deg += log2(3)
        elseif (seq[i] == DNA_K || seq[i] == DNA_Y || seq[i] == DNA_S || seq[i] == DNA_W 
                || seq[i] == DNA_R || seq[i] == DNA_M)
            deg += 1.0
        else # DNA_Gap
            break
        end
    end
    return deg
end

function simulate(prefix, bases, sites; nsims=1000, horizon=length(bases))
    degs = zeros(nsims)
    true_horizon = min(horizon, length(bases))
    for i = 1:nsims
        degs[i] = degeneracy(find_valid_randomer(prefix, bases[1:true_horizon], sites))
    end
    return mean(degs)
end


function *(seq::LongDNASeq, base::DNA)
    seq * LongDNASeq([base])
end


function cutfree_rollout(bases, sites; kwargs...)
    n = length(bases)
    randomer = dna"-" ^ n

    for i = 1:length(sites)
        if ~ispalindromic(sites[i])
            push!(sites, reverse_complement(sites[i]))
        end
    end

    max_len = map(length, sites) |> maximum
    horizon = max_len

    for i = 1:n
        best_base = DNA_Gap
        best_deg = -1.0
        start = max(i - max_len, 1)
        for base in bases[i]
            if ~is_valid(randomer[start:i-1] * base, sites)
                continue
            end
            deg = simulate(randomer[start:i-1] * base, bases[i+1:end], sites; horizon=horizon, kwargs...)
            if deg > best_deg
                best_deg = deg
                best_base = base
            end
        end
        if best_base == DNA_Gap
            break
        else
            randomer[i] = best_base
        end
    end

    return randomer
end

all_bases = LongDNASeq([alphabet(DNA)[2:end]...])
all_ns = [all_bases for _ in 1:10]
all_nas = [dna"NA" for _ in 1:10]
all_N = [dna"N" for _ in 1:10]
sites = [dna"GGTCTC"]

randomer = cutfree_rollout(all_ns, sites, nsims=1000, horizon=6)
deg = degeneracy(randomer)
natdeg = deg / log2(exp(1))
println("Final randomer: $randomer")
println("Final degeneracy: $deg ($natdeg)")

function runtimes(data)
    data[:rollout_objval] = 0.0
    data[:rollout_runtime] = 0.0
    data[:rollout_code] = ""

    data.objval = data.objval ./ log(2)

    bases = [all_bases for _ in 1:20]
    sites = map(x -> split(x, ',') .|> LongDNASeq, data.sites)
    for i = 1:nrow(data)
        stats = @timed cutfree_rollout(bases, sites[i], nsims=1000)
        data[i,:rollout_objval] = degeneracy(stats.value)
        data[i,:rollout_runtime] = stats.time
        data[i,:rollout_code] = convert(String, stats.value)
        println(i)
    end
    return data
end

data = runtimes(CSV.read("./results.csv"))





