function run_Dkl_path_single(i::Int, j::Int, n_sub::Int, seq_start::Array{Int,1}, seq_target::Array{Int,1}, h::Array{Float64,2},J::Array{Float64,4}; q::Int = 21, pc::Float64 = 10^(-5), temp::Float64 = 1.)
    
    L = length(seq_start)
    diff_sites = findall(seq_start .!= seq2) 
    filter!(x -> x !== i, diff_sites)
    filter!(x -> x !== j, diff_sites)
    
    target_i = zeros(q); target_j = zeros(q);
    target_i[seq_target[i]] = 1.; target_j[seq_target[j]] = 1.;
    pseudocount1!(target_i, target_i, pc, q)
    pseudocount1!(target_j, target_j, pc, q)
    
    list_of_seq = [zeros(Int,L) for i in 1:n_sub+1]
    seq = deepcopy(seq_start);
    seq[i] = seq_target[i]; seq[j] = seq_target[j];
    
    for n in 1:n_sub
        costs = []
        for pos in diff_sites
            push!(costs, cost_one_switch(pos, seq, target_i, target_j, i, j, h, J, L; pc = pc, q = q, temp = temp))
        end
        idx = argmin(costs)
        sel_site = diff_sites[idx] 
        seq[sel_site] = seq_target[sel_site]
        list_of_seq[n+1] .= seq 
        filter!(x -> x !== sel_site, diff_sites)
    end
    
    ens = [energy(s, h, J, L) for s in list_of_seq]
    ps_i = [single_site_prob_cond(i, s, h, J, L)[seq_target[i]] for s in list_of_seq]
    ps_j = [single_site_prob_cond(j, s, h, J, L)[seq_target[j]] for s in list_of_seq]
    ps_ij = [dist_proba_2_point(i, j, s, h, J, L, temp)[seq_target[i], seq_target[j]] for s in list_of_seq]
    
    return (seqs = list_of_seq, ens = ens, ps_i = ps_i, ps_j = ps_j, ps_ij = ps_ij)
end



function cost_one_switch(site::Int, seq::Array{Int,1}, seq_target::Array{Int,1},  target_i::Array{Float64,1}, target_j::Array{Float64,1}, i::Int, j::Int, h::Array{Float64,2}, J::Array{Float64,4},L::Int; q::Int = 21, pc::Float64 = 10^(-5), temp::Float64 = 1.)
    
    new_seq = deepcopy(seq)
    new_seq[site] = seq_target[site]
    
    p_i = single_site_prob_cond_with_T(i, new_seq, h, J, L, temp = temp); 
    pseudocount1!(p_i, p_i, pc, q)
    p_j = single_site_prob_cond(j, new_seq, h, J, L);
    pseudocount1!(p_j, p_j, pc, q)
    
    return kldivergence(p_i, target_i) + kldivergence(p_j, target_j)  
end



function pseudocount1!(dest, f1, pc::AbstractFloat, q::Int)
     dest .= ((1-pc) .* f1 ) .+ (pc / q)
end


function run_Dkl_path_double(i::Int, j::Int, n_sub::Int, seq_start::Array{Int,1}, seq_target::Array{Int,1}, h::Array{Float64,2},J::Array{Float64,4}; q::Int = 21, pc::Float64 = 10^(-5), temp::Float64 = 1.)
    
    L = length(seq_start)
    diff_sites = findall(seq_start .!= seq2) 
    filter!(x -> x !== i, diff_sites)
    filter!(x -> x !== j, diff_sites)
    
    target = zeros(q,q); 
    target_[seq_target[i], seq_target[j]] = 1.; 
    pseudocount2!(target, target, pc, q)
    
    list_of_seq = [zeros(Int,L) for i in 1:n_sub+1]
    seq = deepcopy(seq_start);
    seq[i] = seq_target[i]; seq[j] = seq_target[j];
    
    for n in 1:n_sub
        costs = []
        for pos in diff_sites
            push!(costs, cost_double_switch(pos, seq, target, i, j, h, J, L; pc = pc, q = q, temp = temp))
        end
        idx = argmin(costs)
        sel_site = diff_sites[idx] 
        seq[sel_site] = seq_target[sel_site]
        list_of_seq[n+1] .= seq 
        filter!(x -> x !== sel_site, diff_sites)
    end
    
    ens = [energy(s, h, J, L) for s in list_of_seq]
    ps_i = [single_site_prob_cond(i, s, h, J, L)[seq_target[i]] for s in list_of_seq]
    ps_j = [single_site_prob_cond(j, s, h, J, L)[seq_target[j]] for s in list_of_seq]
    ps_ij = [dist_proba_2_point(i, j, s, h, J, L, temp)[seq_target[i], seq_target[j]] for s in list_of_seq]
    
    return (seqs = list_of_seq, ens = ens, ps_i = ps_i, ps_j = ps_j, ps_ij = ps_ij)
end



function cost_double_switch(site::Int, seq::Array{Int,1}, seq_target::Array{Int,1}, target::Array{Float64,2}, i::Int, j::Int, h::Array{Float64,2}, J::Array{Float64,4},L::Int; q::Int = 21, pc::Float64 = 10^(-5), temp::Float64 = 1.)
    
    new_seq = deepcopy(seq)
    new_seq[site] = seq_target[site]
    
    p_ij = dist_proba_2_point(i, j, new_seq, h, J, L, q, temp)     
    pseudocount2!(p_ij, p_ij, pc, q)
    
    return kldivergence(p_ij[:], target[:])  
end

function pseudocount2!(dest, f2, pc::AbstractFloat, q::Int)
     dest .= ((1-pc) .* f2 ) .+ (pc / (q*q))
end

function single_proba_2_point(a::Int, b::Int, i::Int, j::Int, seq::Array{Float64,1}, h::Array{Float64,2}, J::Array{Float64,4}, N::Int, T::Float64)
    log_proba = h[a, i] + h[b, j] + J[a, b, i, j]
    for k in 1:N
        if k !== j
            log_proba += J[a,seq[k], i, k]
        end
    end
    
    for l in 1:N
        if l !== i
            log_proba += J[seq[l], b, l, j]
        end
    end
    return exp(log_proba/T)
end
                            
                            
function dist_proba_2_point(i::Int, j::Int, seq::Array{Float64,1}, h::Array{Float64,2}, J::Array{Float64,4}, N::Int, q::Int, T::Float64) 
    res = zeros(q,q)
    norm = 0.
    for a in 1:q
        for b in 1:q
            res[a,b] = single_proba_2_point(a, b, i, j, seq, h, J, N, T)
            norm += res[a,b]
        end
    end
    return res/norm
end

