function perform_pca(seqs::Array{Int8,2}; q = 21, n_dim = 2)

    W, Meff = compute_weights(seqs, q+1, 0.2)
    L, M = size(seqs)
    
    println("Computing statistics...")
    
    one_hot_seqs = zeros(M,L*(q-1))
    for m in 1:M    # Performs one-hot encoding (ignoring gaps, why?)
        for i in 1:L
            if seqs[i,m] < q
                one_hot_seqs[m,(i-1)*(q-1)+seqs[i,m]] = 1.0
            end
        end
    end

    pij = zeros(L*(q-1),L*(q-1))
    pi = zeros(L*(q-1),)

    for m in 1:M
        for i in 1:L
            if seqs[i,m] < q
                pi[(i-1)*(q-1) + seqs[i,m]] += W[m]
                for j = i:L
                    if seqs[j,m] < q
                        pij[(i-1)*(q-1) + seqs[i,m], (j-1)*(q-1) + seqs[j,m]] += W[m]
                        pij[(j-1)*(q-1) + seqs[j,m], (i-1)*(q-1) + seqs[i,m]] += W[m] 
                    end
                end
            end
        end
    end

    pij /= Meff
    pi /= Meff
    cij = pij - pi*transpose(pi)

    println("Performing eig factorization...")
    F = eigen(cij)
    println("Largest eigenvalues: ", F.values[end:-1:end-n_dim+1])
    #PC = F.vectors[:,end-max_dim+1:end] 

    proj_seqs = zeros(M,n_dim)
    
    for i in 1:n_dim
        proj_seqs[:,i] = one_hot_seqs * F.vectors[:,end-i+1]
    end 
    
    println("Done")
    return proj_seqs,  F.vectors[:, end:-1:end-n_dim+1], F.values[end:-1:end-n_dim+1]/sum(F.values)
end


function get_projection(PC::Array{Float64,2},  seq::Array{Int8,1}; q::Int=21)

    L = length(seq)    
    L_pc, n_dim = size(PC)
	@assert L_pc == L*(q-1)    
    one_hot_sample = zeros(1,L_pc)
    
    # println("Encoding sample...")
    for i in 1:L
        if seq[i] < q # exclude gaps (whatever learning method you used)
            one_hot_sample[1,(i-1)*(q-1)+ seq[i]] = 1.0
        end
    end
	# println("Projecting sample...")
	proj_sample = zeros(n_dim)
    
    for i in 1:n_dim
        proj_sample[i] = one_hot_sample * PC[:,i]
    end 
    
	return proj_sample
end

function get_projection(PC::Array{Float64,2},  seqs::Array{Int8,2}; q::Int=21)

    L, M = size(seqs)    
    L_pc, n_dim = size(PC)
	@assert L_pc == L*(q-1)    
    one_hot_sample = zeros(L_pc, M)
    
    # println("Encoding sample...")
    for m in 1:M
        for i in 1:L
            if seqs[i,m] < q # exclude gaps (whatever learning method you used)
                one_hot_sample[(i-1)*(q-1)+ seqs[i,m], m] = 1.0
            end
        end
    end
	# println("Projecting sample...")
	proj_sample = zeros(M, n_dim)
    
    for i in 1:n_dim
        for m in 1:M
            proj_sample[m,i] = one_hot_sample[:,m]' * PC[:,i]
        end
    end 
    
	return proj_sample
end




function check_pca(folder::String, nat_msa::Array{Int8,2}, sampled_msa::Array{Int8,2})
    pc_nat, PC, expl_var = perform_pca(nat_msa, n_dim = 2);
    pc_sil = get_projection(PC, sampled_msa);
    
    close("all")
    plot_density(pc_nat[:,1], pc_nat[:,2])
    plt.title("Natural sequences")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    savefig(joinpath(folder, "pca_nat.png"))
    
    
    close("all")
    plot_density(pc_sil[:,1], pc_sil[:,2])
    plt.title("Simulated sequences")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    savefig(joinpath(folder, "pca_sim.png"))
    
end
    
    