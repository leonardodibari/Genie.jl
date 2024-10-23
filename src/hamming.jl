function ham_dist(vec1::Array{Int8,1}, vec2::Array{Int8,1})
    return sum(vec1 .!= vec2)
end

function ham_dist(vec::Array{Int8,1}, msa::Array{Int8,2})
    return [ham_dist(vec, msa[:,i]) for i in 1:size(msa,2)]
end

function ham_dist(msa1::Array{Int8,2}, msa2::Array{Int8,2})
    return [ham_dist(msa1[:,i], msa2[:,i]) for i in 1:size(msa1,2)]
end

function ham_dist(vec1::Array{Int,1}, vec2::Array{Int,1})
    return sum(vec1 .!= vec2)
end

function ham_dist(vec::Array{Int,1}, msa::Array{Int,2})
    return [ham_dist(vec, msa[:,i]) for i in 1:size(msa,2)]
end

function ham_dist(msa1::Array{Int,2}, msa2::Array{Int,2})
    return [ham_dist(msa1[:,i], msa2[:,i]) for i in 1:size(msa1,2)]
end

function ham_dist(step_msa::Array{Array{Int8,2},1})
    N_steps = length(step_msa)
    N_chains = size(step_msa[1],2)
    res = zeros(N_steps,N_chains)
    for n in 1:N_steps
        res[n,:] .= ham_dist(step_msa[1], step_msa[n]) 
    end
    return res 
end




function ham_dist_AB_rp2(step_msa::Array{Int8,3}, step_msa_B::Array{Int8,3}, num::Int)
    L, N_points, N_chains = size(step_msa)
    
    res = zeros(N_points)
    res_sq = zeros(N_points)
    for n in 1:N_points
        a = []
        a_sq = []
        for r in 1:num#N_chains
            for p in r+1:num#N_chains
                hh = ham_dist(step_msa[:,n,r], step_msa_B[:,n,p])
                push!(a, hh) 
                push!(a_sq, hh^2)
            end
        end
        res[n] = mean(a)
        res_sq[n] = mean(a_sq)
    end
    return res, res_sq 
end


function ham_dist_AB_rp(step_msa::Array{Array{Int,2},1}, step_msa_B::Array{Array{Int,2},1}, num::Int)
    N_steps = length(step_msa)
    N_chains = size(step_msa[1],2)
    res = zeros(N_steps)
    res_sq = zeros(N_steps)
    for n in 1:N_steps
        a = []
        a_sq = []
        for r in 1:num#N_chains
            for p in r+1:num#N_chains
                hh = ham_dist(step_msa[n][:,r], step_msa_B[n][:,p])
                push!(a, hh) 
                push!(a_sq, hh^2)
            end
        end
        res[n] = mean(a)
        res_sq[n] = mean(a_sq)
    end
    return res, res_sq 
end


function pair_ham_dist_AB_rp(step_msa::Array{Int8,2}, step_msa_B::Array{Int8,2}, num::Int)
    N_chains = size(step_msa,2)
    a = 0; a_sq = 0; tot = num*(num-1)/2
    for r in 1:num#N_chains
        for p in r+1:num#N_chains
            hh = ham_dist(step_msa[:,r], step_msa_B[:,p])
            a +=hh 
            a_sq += hh^2
        end
    end
    return a/tot, a_sq/tot
end
    
function ham_dist_AB_rp3(step_msa::Array{Array{Int8,2},1}, step_msa_B::Array{Array{Int8,2},1}, num::Int)
    N_steps = length(step_msa)
    N_chains = size(step_msa[1],2)
    res = zeros(N_steps)
    res_sq = zeros(N_steps)
    for n in 1:N_steps
        a,a_sq = pair_ham_dist_AB_rp(step_msa[n], step_msa_B[n],num)
        res[n] = a
        res_sq[n] = a_sq
    end
    return res, res_sq 
end

function ham_dist_AB_rp(step_msa::Array{Array{Int8,2},1}, step_msa_B::Array{Array{Int8,2},1}, num::Int)
    N_steps = length(step_msa)
    N_chains = size(step_msa[1],2)
    res = zeros(N_steps)
    res_sq = zeros(N_steps)
    for n in 1:N_steps
        a = []
        a_sq = []
        for r in 1:num#N_chains
            for p in r+1:num#N_chains
                hh = ham_dist(step_msa[n][:,r], step_msa_B[n][:,p])
                push!(a, hh) 
                push!(a_sq, hh^2)
            end
        end
        res[n] = mean(a)
        res_sq[n] = mean(a_sq)
    end
    return res, res_sq 
end


function ham_dist(step_msa::Array{Array{Int,2},1})
    N_steps = length(step_msa)
    N_chains = size(step_msa[1],2)
    res = zeros(N_steps,N_chains)
    for n in 1:N_steps
        res[n,:] .= ham_dist(step_msa[1], step_msa[n]) 
    end
    return res 
end

function pairwise_ham_dist(msa::Array{Int8,2}; n_seq = 100, all = false)
    res = []
    for i in 1:n_seq
        for j in i+1:n_seq
            push!(res,ham_dist(msa[:,i], msa[:,j]))
        end
    end
    if all == false
        return mean(res)
    else 
        return res
    end
end