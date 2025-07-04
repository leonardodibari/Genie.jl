
function insert_muts(muts, seq_start, seq_target, i, j, h, J; q = 21)
    L = length(seq_start);
    seq = deepcopy(seq_start)
    seq[i] = seq_target[i]
    seq[j] = seq_target[j]
    ps = [];
    for m in muts
        seq[m] = seq_target[m]
        push!(ps, dist_proba_2_point(i, j, seq, h, J, L, q, 1.)[seq_target[i], seq_target[j]])
    end
        
    return ps
end
    
  
function try_order_muts(muts, n_shuffles, seq_start, seq_target, i, j, h, J; q = 21)
    res = []
    pps = []
    for n in 1:n_shuffles
        new_muts = shuffle(muts)
        push!(res, new_muts)
        push!(pps, insert_muts(new_muts, seq_start, seq_target, i, j, h, J; q = 21))
    end
        
    return hcat(res...), hcat(pps...)
end
        
    

    


function run_Dkl_path_single(i::Int, n_sub::Int, seq_start::Array{Int,1}, seq_target::Array{Int,1}, h::Array{Float64,2},J::Array{Float64,4}; q::Int = 21, pc::Float64 = 10^(-5), temp::Float64 = 1., reg = 0., reg_type = "dE")
    
    L = length(seq_start)
    diff_sites = findall(seq_start .!= seq_target) 
    filter!(x -> x !== i, diff_sites)
    #filter!(x -> seq_start[x] !== 21, diff_sites)
    #filter!(x -> seq_target[x] !== 21, diff_sites)
    
    println(length(diff_sites))
    
    
    target_i = zeros(q); 
    target_i[seq_target[i]] = 1.; 
    pseudocount1!(target_i, target_i, pc, q)
    
    list_of_seq = [zeros(Int,L) for i in 1:n_sub+1]
    list_of_seq[1] .= seq_start
    seq = deepcopy(seq_start);
    #seq[i] = seq_target[i]; 
    
    for n in 1:n_sub
        costs = zeros(length(diff_sites))
        costs_energy = zeros(length(diff_sites))
        @tasks for ii in 1:length(diff_sites)
            pos = diff_sites[ii]
            costs[ii] =  cost_one_switch(pos, seq, seq_target, target_i, i, h, J, L; pc = pc, q = q, temp = temp)
            costs_energy[ii] = single_mut_dE(seq, h, J, seq_target[pos], pos, L)
        end
        
        #=println(n)
        println(size(costs_energy))
        println(size((reg .* relu.(costs_energy))))
        println(size(costs))=#
        
        if reg_type == "dE"
            idx = argmin(costs .+ (reg .* costs_energy))
        elseif reg_type == "abs"
            idx = argmin(costs .+ (reg .* abs.(costs_energy)))
        elseif reg_type == "relu"
            idx = argmin(costs .+ (reg .* relu.(costs_energy)))
        end
                       
        sel_site = diff_sites[idx] 
        seq[sel_site] = seq_target[sel_site]
        list_of_seq[n+1] .= seq 
        filter!(x -> x !== sel_site, diff_sites)
        
        ppp = single_site_prob_cond_with_T(i, seq, h, J, L, temp = temp)
        if ppp[seq_target[i]] > ppp[seq_start[i]]
            seq[i] = seq_target[i]
        end 
    end
    
    ens = [energy(s, h, J, L) for s in list_of_seq]
    ps_i = [single_site_prob_cond(i, s, h, J, L)[seq_target[i]] for s in list_of_seq]
    
    return (seqs = list_of_seq, ens = ens, ps_i = ps_i)
end



function cost_one_switch(site::Int, seq::Array{Int,1}, seq_target::Array{Int,1},  target_i::Array{Float64,1}, i::Int, h::Array{Float64,2}, J::Array{Float64,4},L::Int; q::Int = 21, pc::Float64 = 10^(-5), temp::Float64 = 1.)
    
    new_seq = deepcopy(seq)
    new_seq[site] = seq_target[site]
    
    p_i = single_site_prob_cond_with_T(i, new_seq, h, J, L, temp = temp); 
    pseudocount1!(p_i, p_i, pc, q)
    
    
    return kldivergence(p_i, target_i)
end



function pseudocount1!(dest, f1, pc::AbstractFloat, q::Int)
     dest .= ((1-pc) .* f1 ) .+ (pc / q)
end


function run_Dkl_path_double(i::Int, j::Int, n_sub::Int, seq_start::Array{Int,1}, seq_target::Array{Int,1}, h::Array{Float64,2},J::Array{Float64,4}; q::Int = 21, pc::Float64 = 10^(-5), temp::Float64 = 1.)
    
    L = length(seq_start)
    diff_sites = findall(seq_start .!= seq_target) 
    filter!(x -> x !== i, diff_sites)
    filter!(x -> x !== j, diff_sites)
    filter!(x -> seq_start[x] !== 21, diff_sites)
    filter!(x -> seq_target[x] !== 21, diff_sites)
    
    
    println(length(diff_sites))
    
    target = zeros(q,q); 
    target[seq_target[i], seq_target[j]] = 1.; 
    pseudocount2!(target, target, pc, q)
    
    list_of_seq = [zeros(Int,L) for i in 1:n_sub+1]
    list_of_seq[1] .= seq_start
    seq = deepcopy(seq_start);
    seq[i] = seq_target[i]; seq[j] = seq_target[j];
    
    muts = []
    
    for n in 1:n_sub
        costs = []
        for pos in diff_sites
            push!(costs, cost_double_switch(pos, seq, seq_target, target, i, j, h, J, L; pc = pc, q = q, temp = temp))
        end
        idx = argmin(costs)
        sel_site = diff_sites[idx] 
        seq[sel_site] = seq_target[sel_site]
        list_of_seq[n+1] .= seq 
        filter!(x -> x !== sel_site, diff_sites)
        push!(muts, sel_site)
    end
    
    ens = [energy(s, h, J, L) for s in list_of_seq]
    ps_i = [single_site_prob_cond(i, s, h, J, L)[seq_target[i]] for s in list_of_seq]
    ps_j = [single_site_prob_cond(j, s, h, J, L)[seq_target[j]] for s in list_of_seq]
    ps_ij = [dist_proba_2_point(i, j, s, h, J, L, q, 1.)[seq_target[i], seq_target[j]] for s in list_of_seq]
    
    return (seqs = list_of_seq, ens = ens, ps_i = ps_i, ps_j = ps_j, ps_ij = ps_ij, muts = muts)
end


function run_pij_path_double(i::Int, j::Int, n_sub::Int, seq_start::Array{Int,1}, seq_target::Array{Int,1}, h::Array{Float64,2},J::Array{Float64,4}; q::Int = 21, pc::Float64 = 10^(-5), temp::Float64 = 1.)
    
    L = length(seq_start)
    diff_sites = findall(seq_start .!= seq_target) 
    filter!(x -> x !== i, diff_sites)
    filter!(x -> x !== j, diff_sites)
    filter!(x -> seq_start[x] !== 21, diff_sites)
    filter!(x -> seq_target[x] !== 21, diff_sites)
    
    
    println(length(diff_sites))
    
    target = zeros(q,q); 
    target[seq_target[i], seq_target[j]] = 1.; 
    pseudocount2!(target, target, pc, q)
    
    list_of_seq = [zeros(Int,L) for i in 1:n_sub+1]
    list_of_seq[1] .= seq_start
    seq = deepcopy(seq_start);
    seq[i] = seq_target[i]; seq[j] = seq_target[j];
    
    muts = []
    
    for n in 1:n_sub
        costs = []
        for pos in diff_sites
            new_seq = deepcopy(seq)
            new_seq[pos] = seq_target[pos]
            push!(costs, dist_proba_2_point(i, j, new_seq, h, J, 
                    L, q, temp)[seq_target[i],seq_target[j]])
        end
        idx = argmax(costs)
        sel_site = diff_sites[idx] 
        seq[sel_site] = seq_target[sel_site]
        list_of_seq[n+1] .= seq 
        filter!(x -> x !== sel_site, diff_sites)
        push!(muts, sel_site)
    end
    
    ens = [energy(s, h, J, L) for s in list_of_seq]
    ps_i = [single_site_prob_cond(i, s, h, J, L)[seq_target[i]] for s in list_of_seq]
    ps_j = [single_site_prob_cond(j, s, h, J, L)[seq_target[j]] for s in list_of_seq]
    ps_ij = [dist_proba_2_point(i, j, s, h, J, L, q, 1.)[seq_target[i], seq_target[j]] for s in list_of_seq]
    
    return (seqs = list_of_seq, ens = ens, ps_i = ps_i, ps_j = ps_j, ps_ij = ps_ij, muts = muts)
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

function single_proba_2_point(a::Int, b::Int, i::Int, j::Int, seq::Array{Int,1}, h::Array{Float64,2}, J::Array{Float64,4}, N::Int, T::Float64)
    log_proba = h[a, i] + h[b, j] + J[a, i, b, j]
    for k in 1:N
        if k !== j
            log_proba += J[a, i, seq[k], k]
        end
    end
    
    for l in 1:N
        if l !== i
            log_proba += J[seq[l], l, b, j]
        end
    end
    return exp(log_proba/T)
end
                            
                            
function dist_proba_2_point(i::Int, j::Int, seq::Array{Int,1}, h::Array{Float64,2}, J::Array{Float64,4}, N::Int, q::Int, T::Float64) 
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







function run_Dkl_path_energy(i::Int, j::Int, n_sub::Int, seq_start::Array{Int,1}, seq_target::Array{Int,1}, h::Array{Float64,2},J::Array{Float64,4}; q::Int = 21, pc::Float64 = 10^(-5), temp::Float64 = 1.)
    
    L = length(seq_start)
    diff_sites = findall(seq_start .!= seq_target) 
    filter!(x -> x !== i, diff_sites)
    filter!(x -> x !== j, diff_sites)
    filter!(x -> seq_start[x] !== 21, diff_sites)
    filter!(x -> seq_target[x] !== 21, diff_sites)
    
    target = zeros(q,q); 
    target[seq_target[i], seq_target[j]] = 1.; 
    pseudocount2!(target, target, pc, q)
    
    list_of_seq = [zeros(Int,L) for i in 1:n_sub+1]
    list_of_seq[1] .= seq_start
    seq = deepcopy(seq_start);
    seq[i] = seq_target[i]; seq[j] = seq_target[j];
    
    for n in 1:n_sub
        costs = []
        for pos in diff_sites
            push!(costs, cost_energy_switch(pos, seq, seq_target, i, j, h, J, L; pc = pc, q = q, temp = temp))
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
    ps_ij = [dist_proba_2_point(i, j, s, h, J, L, q, 1.)[seq_target[i], seq_target[j]] for s in list_of_seq]
    
    return (seqs = list_of_seq, ens = ens, ps_i = ps_i, ps_j = ps_j, ps_ij = ps_ij)
end




function cost_energy_switch(site::Int, seq::Array{Int,1}, seq_target::Array{Int,1}, i::Int, j::Int, h::Array{Float64,2}, J::Array{Float64,4},L::Int; q::Int = 21, pc::Float64 = 10^(-5), temp::Float64 = 1.)
            
    return single_mut_dE(seq, h, J, seq_target[site], site, L)  
end



function run_Dkl_path_mixed(i::Int, j::Int, n_sub::Int, seq_start::Array{Int,1}, seq_target::Array{Int,1}, h::Array{Float64,2},J::Array{Float64,4}; q::Int = 21, pc::Float64 = 10^(-5), temp::Float64 = 1., reg = 10^-2)
    
    L = length(seq_start)
    diff_sites = findall(seq_start .!= seq_target) 
    filter!(x -> x !== i, diff_sites)
    filter!(x -> x !== j, diff_sites)
    filter!(x -> seq_start[x] !== 21, diff_sites)
    filter!(x -> seq_target[x] !== 21, diff_sites)
    
    target = zeros(q,q); 
    target[seq_target[i], seq_target[j]] = 1.; 
    pseudocount2!(target, target, pc, q)
    
    list_of_seq = [zeros(Int,L) for i in 1:n_sub+1]
    list_of_seq[1] .= seq_start
    seq = deepcopy(seq_start);
    seq[i] = seq_target[i]; seq[j] = seq_target[j];
    
    muts = []
    
    for n in 1:n_sub
        cost_energy = []
        cost_p = []
        for pos in diff_sites 
            push!(cost_energy, cost_energy_switch(pos, seq, seq_target, i, j, h, J, L; pc = pc, q = q, temp = temp))
            push!(cost_p, cost_double_switch(pos, seq, seq_target, target, i, j, h, J, L; pc = pc, q = q, temp = temp))
        end
        idx = argmin(cost_p .+  (reg .* cost_energy) )
        sel_site = diff_sites[idx] 
        seq[sel_site] = seq_target[sel_site]
        list_of_seq[n+1] .= seq 
        filter!(x -> x !== sel_site, diff_sites)
        push!(muts, sel_site)
    end
    
    ens = [energy(s, h, J, L) for s in list_of_seq]
    ps_i = [single_site_prob_cond(i, s, h, J, L)[seq_target[i]] for s in list_of_seq]
    ps_j = [single_site_prob_cond(j, s, h, J, L)[seq_target[j]] for s in list_of_seq]
    ps_ij = [dist_proba_2_point(i, j, s, h, J, L, q, 1.)[seq_target[i], seq_target[j]] for s in list_of_seq]
    
    return (seqs = list_of_seq, ens = ens, ps_i = ps_i, ps_j = ps_j, ps_ij = ps_ij, muts = muts)
end

