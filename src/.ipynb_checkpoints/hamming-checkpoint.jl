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

function ham_dist(step_msa::Array{Array{Int,2},1})
    N_steps = length(step_msa)
    N_chains = size(step_msa[1],2)
    res = zeros(N_steps,N_chains)
    for n in 1:N_steps
        res[n,:] .= ham_dist(step_msa[1], step_msa[n]) 
    end
    return res 
end

function pairwise_ham_dist(msa::Array{Int8,2}; n_seq = 100)
    res = []
    msa = 
    for i in 1:n_seq
        for j in i+1:n_seq
            push!(res,ham_dist(msa[:,i], msa[:,j]))
        end
    end
    return mean(res)
end