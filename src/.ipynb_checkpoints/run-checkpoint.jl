
# run.jl
# Main script to execute and manage the analyses for energy and timescale calculations.

"""
    run_evolution(start_msa, h::Array{T,2}, J::Array{T,4};
                   N_steps::Int = 100, 
                   temp::Float64 = 1.0,  
                   p::Float64 = 0.5, 
                   N_points::Union{Int, Nothing} = nothing, 
                   each_step::Union{Int, Nothing} = nothing, 
                   rand_init::Bool = false, 
                   q::Int = 21, 
                   codon_bias::Union{Nothing, Dict{String, Float64}} = nothing, 
                   verbose::Bool = false) where {T}

Run the evolutionary simulation of multiple sequence alignments (MSA) using Gibbs sampling and Metropolis algorithms.

# Arguments
- `start_msa`: Initial multiple sequence alignment (MSA) as a matrix. *(non-optional)*
- `h`: A 2D array of external field parameters. *(non-optional)*
- `J`: A 4D array of interaction parameters, should be of shape `(q, L, q, L)`. *(non-optional)*
- `N_steps`: Total number of evolutionary steps to simulate (default is `100`). *(optional)*
- `temp`: Temperature parameter for sampling algorithms (default is `1.0`). *(optional)*
- `p`: Probability of using Gibbs sampling over Metropolis sampling (default is `0.5`). *(optional)*
- `N_points`: Number of points to collect along the trajectory in logarithmic scale (default is `nothing`). *(optional)*
- `each_step`: Step interval for collecting MSAs (default is `nothing`). *(optional)*
- `rand_init`: Initialize chains randomly if set to `true` (default is `false`). *(optional)*
- `q`: Number of different codons (default is `21`). *(optional)*
- `codon_bias`: Dictionary specifying codon usage bias (default is `nothing`). *(optional)*
- `verbose`: Print detailed output during execution if set to `true` (default is `false`). *(optional)*

# Returns
- A tuple containing:
  - `step_msa`: List of MSAs collected at specified steps (if `N_points` or `each_step` is specified).
  - `step_msa_dna`: Corresponding DNA sequences.
  - `codon_usage`: Dictionary of codon usage.
  - `steps`: The list of steps at which MSAs were collected.
  - `msa`: Final MSA collected after all steps (if neither `N_points` nor `each_step` is specified).
  - `msa_dna`: Final DNA sequences corresponding to the MSA.

# Errors
- Throws errors for invalid dimensions of input arrays `J` and `h`, conflicts between `N_points` and `each_step`, only one of the two must be different from nothing and they must be smaller than `N_steps`.
"""


function run_evolution(start_msa, h::Array{T,2}, J::Array{T,4}; 
        N_steps::Int = 100, 
        temp = 1.0,  
        p = 0.5, 
        N_points::Union{Int, Nothing} = nothing, 
        each_step::Union{Int, Nothing} = nothing, 
        rand_init = false, 
        q = 21, 
        codon_bias::Union{Nothing, Dict{String, Float64}} = nothing, 
        verbose = false) where {T}
    
    
    L, N_chains = size(start_msa)
    if (size(J,1) !== size(J,3)) || (size(J,2) !== size(J,4))
        error("Size of J should be (q,L,q,L)")
    elseif (size(J,2) !== L) || (size(h,2) !== L)
        error("Length of sequences different from length of parameters")
    elseif N_points !== nothing && each_step !== nothing
        error("You cannot specify both N_points and each_step: 
            \n-if you specify N_points, you get some msa (N_points) in logarithmic scale along the trajectory
            \n-if you specify each_step, you get a msa every each_step along the trajectory
            \n-if you don't specify anything then you get only the msa at the final step \n")
    end
    
    rng = random_gens(N_chains)
    codon_net = create_nested_codon_dict()
    length_of_moves = create_length_dict(codon_net)
    count = 0
    temp = T(temp)
    all_codons = vcat([amino2cod[i] for i in 1:20]...)
    push!(all_codons, "TAG")
    push!(all_codons, "TAA")
    push!(all_codons, "TGA")
    #prob_all_cod = [1/64 for _ in 1:length(all_codons)]
    
    if codon_bias == nothing
        no_cod_bias = Dict(x => T(1/length(amino2cod[cod2amino[x]])) for x in keys(cod2amino))
        codon_usage = no_cod_bias
    else
        codon_usage = codon_bias
    end
    
    if rand_init == true
        println("Random Initialization")
        chains = [Chain(Int8.(rand(1:q, L)), q, rng[n]) for n in 1:N_chains]
    else
        chains = [Chain(start_msa[:,n], q, rng[n]) for n in 1:N_chains]
    end
     
    if N_points !== nothing 
        if N_points > N_steps
            error("N_points must be smaller than N_steps")
        end
        steps = unique([trunc(Int,10^y) for y in range(log10(1), log10(N_steps), length=N_points)])
        step_msa = [zeros(Int8, (L, N_chains)) for i in 1:length(steps)]
        step_msa_dna = [Matrix{String}(undef, L, N_chains) for i in 1:length(steps)]
    end
    
    if each_step !== nothing 
        if each_step > N_steps
            error("each_step must be smaller than N_steps")
        end
        steps = [i for i in 1:each_step:N_steps]
        step_msa = [zeros(Int8, (L, N_chains)) for i in 1:length(steps)]
        step_msa_dna = [Matrix{String}(undef, L, N_chains) for i in 1:length(steps)]
    end
    
     
    @inbounds for t in 1:N_steps
        
        if ((N_points !== nothing) || (each_step !== nothing)) && (t in steps)
            count += 1
            if verbose == true
                println(t)
            end
            @tasks for n in 1:N_chains
                for i in 1:L
                    step_msa[count][i,n] = chains[n].seq[i]
                    step_msa_dna[count][i,n] = chains[n].DNA[i]
                end
            end
        end
        
        if rand() < 1 - p
            #sampling gibbs with probability p
            run_gibbs_sampling!(chains, h, J, codon_net, codon_usage, length_of_moves, N_chains, temp, L)
        else
            #sampling metropolis with probability 1-p
            run_metropolis_indels!(chains, h, J, all_codons, codon_usage, N_chains, temp,L)
        end 
            
    end  
    
    if (N_points !== nothing) || (each_step !== nothing) 
        return (step_msa = step_msa, step_msa_dna = step_msa_dna, codon_usage = codon_usage, steps = steps, p = p, temp = temp)
    else
        msa = Int8.(zeros(L, N_chains))
        msa_dna = Matrix{String}(undef, L, N_chains)
        for n in 1:N_chains
            for i in 1:L
                msa[i,n] = Int8.(chains[n].seq[i])
                msa_dna[i,n] = chains[n].DNA[i]
            end
        end
        return (msa = msa, msa_dna = msa_dna, codon_usage = codon_usage, p = p, temp = temp)
    end
    
end
    


