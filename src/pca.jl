

function one_hot(array::Array{Int8, 2})
    n = size(array, 1)
    m = size(array, 2)
    one_hot_array = zeros(Int8, n, m * 20)
    for i in 1:n
        for j in 1:m
            if array[i,j] == 21
                # Gap amino acid
                index = (j-1) * 20 + 1
                one_hot_array[i,index:index+19] = zeros(Int8, 20)
            else
                # Non-gap amino acid
                index = (j-1) * 20 + array[i,j]
                one_hot_array[i,index] = 1
            end
        end
    end
    return one_hot_array
end


function sample_and_map(arnet, M; nsample = 10000)
    Zs_one_hot = one_hot(sample(arnet, nsample))

    return predict(M, Zs_one_hot);
end

Z,_ = AttentionDCA.quickread(fastas[i])
Z_one_hot = one_hot(Z)
M = fit(PCA, Z_one_hot, maxoutdim=2);

Z_pc = predict(M, Z_one_hot);

#Sample from an AttentionDCA
arnet = arouts[i][1]
Zs_aratt_pc = sample_and_map(arnet, M, nsample = 2*size(Z_one_hot,2));


j = 3
lim = 5.9
fig, ax = subplots(1, figsize = (5,5))
#rp = rand(1:size(Z,2),50000)
#plot_density(Zs[j][1,rp], Zs[j][2,rp])
plot_density(Zs[j][1,:], Zs[j][2,:])
xticks(fontsize = 15)
yticks(fontsize = 15)
ax.set_xlim([-lim,lim])
ax.set_ylim([-lim,lim])
ax.set_xlabel("First principal component", fontsize = 15)
ax.set_ylabel("Second principal component", fontsize = 15)
ax.set_title(names[i]*", "*versions[j], fontsize = 15)
ax.set_aspect("equal")
savefig("../graphs/"names[i]"pca"vers[j]".pdf", bbox_inches="tight")
