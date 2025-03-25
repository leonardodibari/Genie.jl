

function prob_mut_one_site(mu::T; q::Int = 21) where {T}
    return T(((q-1)/q)*(1-exp(-mu)))   #jukes cantor model
end

   
    
function random_mutations!(msa::Array{Int8,2}, abund_dict::Dict{Vector{Int8}, Int64}, mu::T; q::Int = 21, verbose = false) where {T}
    
    L,M = size(msa)
    distrib = Binomial(L, prob_mut_one_site(mu, q = q))
    
    for m in 1:M
        K = rand(distrib) #number of mutations in sequence 
        aa = rand(1:q, K)
        i = sample(1:L, K, replace = false)
        for k in 1:K
            msa[i[k],m] = aa[k]
        end
    end
    
    empty!(abund_dict)
    seq_count_dict!(abund_dict, msa)
    #=cc = 0
    for val in values(abund_dict)
        cc += val
    end
    println(cc)=#
    clean_dict!(abund_dict)
    #=cc = 0
    for val in values(abund_dict)
        cc += val
    end
    println(cc)=#
      
    if verbose == true
        println("$(100*length(abund_dict)/M) % unique seqs after mutagenesis")
    end
end
    

function selection!(h::Array{T,2}, J::Array{T,4}, abund_dict::Dict{Vector{Int8}, Int64}, L::Int, M::Int; temp::T = 1., mu_bind::T = 18.6, verbose = false) where {T}
        
    surv_prob = []
    for seq in keys(abund_dict)
        a = exp( - ((energy(seq, h, J, L) - mu_bind ) / temp))
        push!(surv_prob, a/(1+a))    
    end
        
    #println("Extrema of surv_prob is $(extrema(surv_prob))")
    
    idx = 1; 
    for seq in keys(abund_dict)
        distrib = Binomial(abund_dict[seq], surv_prob[idx])
        abund_dict[seq] = rand(distrib)
        idx+=1    
    end
    
    clean_dict!(abund_dict)
    
    #=cc = 0
    for val in values(abund_dict)
        cc += val
    end
    println(cc)=#
    if verbose == true    
        println("$(100*length(abund_dict)/M) % unique seqs after selection")
    end
end


function amplification!(final_msa::Array{Int8,2}, abund_dict::Dict{Vector{Int8}, Int64})
    
    L, Mf = size(final_msa);
        
    seqs = [seq for seq in keys(abund_dict)];   
    w = [abund_dict[seq] for seq in keys(abund_dict)];
    #w2 = [abund_dict[seq] for seq in seqs];
    
    M = length(seqs);
    
    #println(sum(abs2, w .- w2))
    
    idxs = sample(1:M, Weights(w), Mf);
    
    for m in 1:Mf
        for i in 1:L
            final_msa[i,m] = seqs[idxs[m]][i]
        end
    end
    
end


function seq_count_dict!(dict::Dict{Vector{Int8}, Int64}, msa::Array{Int8,2})

    M = size(msa,2)
    
    for m=1:M
        haskey(dict,msa[:,m]) ?  dict[msa[:,m]]+=1 : dict[msa[:,m]]=1
    end

    return dict

end



function clean_dict!(d::Dict{Vector{Int8}, Int64})
    for k in keys(d)  # Iterate over keys
        if d[k] == 0
            delete!(d, k)  # Remove key with value 0
        end
    end
end
    

function run_dir_evol(start_msa::Array{Int8,2}, h::Array{T,2}, J::Array{T,4};
                   rounds::Int = 4, 
                   each_step::Union{Int, Nothing} = nothing,
                   seq_steps::Union{Int, Nothing} = nothing,
                   seq_reads::Int = 100,
                   temp::Float64 = 1.0,  
                   mu::T = 10^-2,
                   mu_bind::Float64 = 18.6,
                   q::Int = 21,
                   verbose = false) where {T}
    
    L,M = size(start_msa);
    
    if M < seq_reads
        error("Increase number of sequences, seq_reads must be smaller than amount of sequeces")
    end
    
    if (size(J,1) !== size(J,3)) || (size(J,2) !== size(J,4))
        error("Size of J should be (q,L,q,L)")
    elseif (size(J,2) !== L) || (size(h,2) !== L)
        error("Length of sequences different from length of parameters")
    elseif (size(J,1) !== size(h,1)) || (size(J,3) !== size(h,1))
        error("Issues with color size")
    elseif seq_steps !== nothing && each_step !== nothing
        error("You cannot specify both seq_steps and each_step: 
            \n-if you specify seq_steps, you get a msa on the times you suggested along the trajectory
            \n-if you specify each_step, you get a msa every each_step along the trajectory
            \n-if you don't specify anything then you get only the msa at the final round \n")
    end
        
    
    if each_step !== nothing
        if each_step > rounds
            error("each_step must be smaller than rounds")
        end
        steps = [i for i in 1:each_step:rounds] 
        step_msa = [zeros(Int8, (L, seq_reads)) for i in 1:length(steps)]
    end
    
    if seq_steps !== nothing
        if maximum(seq_steps) > rounds
            error("Sequencing steps must be smaller than rounds")
        elseif is_positive_increasing(seq_steps) == false
            error("Sequencing steps must be strictly positive and increasing")
        end
        steps = deepcopy(seq_steps)
        step_msa = [zeros(Int8, (L, seq_reads)) for i in 1:length(steps)]
    end
    
    
    final_msa = deepcopy(start_msa);
    abund_dict = Dict{Vector{Int8},Int}()
    count = 0
    for r in 1:rounds
        if verbose == true
            println()
            println("Round $(r)")
        end
        Genie.random_mutations!(final_msa, abund_dict, mu, verbose = verbose);
        Genie.selection!(h, J, abund_dict, L, M, temp = temp, mu_bind = mu_bind, verbose = verbose);
        clean_dict!(abund_dict)
        Genie.amplification!(final_msa, abund_dict) 
        
        if ((seq_steps !== nothing) || (each_step !== nothing)) && (r in steps)
            count+=1
            step_msa[count] .= final_msa[:,rand(1:M,seq_reads)]
        end
        
    end
    
    if (seq_steps !== nothing) || (each_step !== nothing)
        return (start_msa = start_msa, steps = steps, step_msa = step_msa, mu = mu, mu_bind = mu_bind, temp = temp)
    else
        return (start_msa = start_msa, rounds = rounds, final_msa = final_msa[:,rand(1:M,seq_reads)], mu = mu, mu_bind = mu_bind, temp = temp)
    end
end


function run_dir_evol(start_seq, N_start::Int, h::Array{T,2}, J::Array{T,4};
                   rounds::Int = 4, 
                   each_step::Union{Int, Nothing} = nothing,
                   seq_steps::Union{Int, Nothing} = nothing,
                   seq_reads::Int = 100,
                   temp::Float64 = 1.0,  
                   mu::T = 10^-2,
                   mu_bind::Float64 = 18.6,
                   q::Int = 21,
                   verbose = false) where {T}
    
    
    start_msa = Int8.(hcat([start_seq for i in 1:N_start]...));
    return run_dir_evol(start_msa, h, J,
                   rounds = rounds,
                   each_step = each_step,
                   seq_reads = seq_reads,
                   temp = temp,  
                   mu = mu,
                   mu_bind = mu_bind,
                   q = q,
                   verbose = verbose) 
end


