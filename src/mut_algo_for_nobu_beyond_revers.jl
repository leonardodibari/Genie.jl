

function run_Dkl_path_double_beyond_rev(i::Int, j::Int, n_sub::Int, seq_start::Array{Int,1}, seq_target::Array{Int,1}, h::Array{Float64,2},J::Array{Float64,4}; q::Int = 21, pc::Float64 = 10^(-5), temp::Float64 = 1., reg = 0.)
    
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
        costs = zeros(21,length(diff_sites))
        cost_energy = zeros(21,length(diff_sites))
        @tasks for pos in 1:length(diff_sites)
            for aa in 1:21
                if aa !== seq[diff_sites[pos]]
                    costs[aa,pos] = cost_double_switch_beyond_rev(diff_sites[pos], seq, aa, target, i, j, h, J, L; pc = pc, q = q, temp = temp)
                    cost_energy[aa,pos] = cost_energy_switch_beyond_rev(pos, seq, aa, i, j, h, J, L; pc = pc, q = q, temp = temp)
                end
            end
            costs[seq[diff_sites[pos]], pos] = 10000.
            cost_energy[seq[diff_sites[pos]], pos] = 10000.
        end
        aa,idx = Tuple(argmin(costs .+ (reg .* cost_energy) ))
        sel_site = diff_sites[idx] 
        seq[sel_site] = aa
        list_of_seq[n+1] .= seq 
        filter!(x -> x !== sel_site, diff_sites)
    end
    
    ens = [energy(s, h, J, L) for s in list_of_seq]
    ps_i = [single_site_prob_cond(i, s, h, J, L)[seq_target[i]] for s in list_of_seq]
    ps_j = [single_site_prob_cond(j, s, h, J, L)[seq_target[j]] for s in list_of_seq]
    ps_ij = [dist_proba_2_point(i, j, s, h, J, L, q, 1.)[seq_target[i], seq_target[j]] for s in list_of_seq]
    
    return (seqs = list_of_seq, ens = ens, ps_i = ps_i, ps_j = ps_j, ps_ij = ps_ij)
end



function cost_double_switch_beyond_rev(site::Int, seq::Array{Int,1}, aa::Int, target::Array{Float64,2}, i::Int, j::Int, h::Array{Float64,2}, J::Array{Float64,4},L::Int; q::Int = 21, pc::Float64 = 10^(-5), temp::Float64 = 1.)
    
    new_seq = deepcopy(seq)
    new_seq[site] = aa
    
    p_ij = dist_proba_2_point(i, j, new_seq, h, J, L, q, temp)     
    pseudocount2!(p_ij, p_ij, pc, q)
    
    return kldivergence(p_ij[:], target[:])  
end



function cost_energy_switch_beyond_rev(site::Int, seq::Array{Int,1},aa::Int, i::Int, j::Int, h::Array{Float64,2}, J::Array{Float64,4},L::Int; q::Int = 21, pc::Float64 = 10^(-5), temp::Float64 = 1.)
            
    return single_mut_dE(seq, h, J, aa, site, L)  
end

