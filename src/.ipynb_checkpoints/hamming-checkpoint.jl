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