function G4(msa::Array{Int8,2}, start_wt::Array{Int8,1})
    
    L, N_chains = size(msa)
    
    similarity = (msa .== start_wt)
    single_mean = mean(similarity, dims = 2)[:]
    res = zeros(L,L)
    
    for i in 1:L
        res[i,i] = single_mean[i]*(1-single_mean[i])
        for j in i+1:L
            res[i,j] = mean(similarity[i,:] .* similarity[j,:]) - (single_mean[i]*single_mean[j])          
            res[j,i] = res[i,j]
        end
    end
    
    return res
    
end


function G4(msa::Array{Int8,2}, start_wt::Array{Int8,1}, w)
    
    L, N_chains = size(msa)
    
    similarity = (msa .== start_wt)
    single_mean = [mean(similarity[i,:], Weights(w)) for i in 1:L]
    res = zeros(L,L)
    
    for i in 1:L
        res[i,i] = single_mean[i]*(1-single_mean[i])
        for j in i+1:L
            res[i,j] = mean(similarity[i,:] .* similarity[j,:], Weights(w)) - (single_mean[i]*single_mean[j])         
            res[j,i] = res[i,j]
        end
    end
    
    
    
    return res
    
end

            
            