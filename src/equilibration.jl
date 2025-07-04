function correlations(nf;
    lambdaJ::Real=0.0001,
    lambdaH::Real=0.0001,
    theta::Union{Symbol,Real}=:auto,
    max_gap_fraction::Real=0.9,
    remove_dups::Bool=true,
    q=21,
    nsamples=10_000)

    Wnat, Znat, _,_,_ = ArDCA.read_fasta(nf, max_gap_fraction, theta, remove_dups)
    L = size(Znat,1)
    #println(L)
    arnet, arvar=ardca(nf, lambdaJ=lambdaJ, lambdaH=lambdaH)
    Zgen = ArDCA.sample(arnet, nsamples)
    Pinat, Pijnat = DCAUtils.compute_weighted_frequencies(Int8.(Znat), Wnat, q + 1)
    Pigen, Pijgen = DCAUtils.compute_weighted_frequencies(Int8.(Zgen), 1.0 ./ ones(nsamples), q + 1)
    Cgen  = reshape(Pijgen - Pigen * Pigen',q,L,q,L)
    Cnat  = reshape(Pijnat - Pinat * Pinat',q,L,q,L)
    Cnatod = Float64[]
    Cgenod = Float64[]
    for a =1:q
        for i=1:L-1
            for b=1:q
                for j=i+1:L
                    push!(Cnatod, Cnat[a,i,b,j])
                    push!(Cgenod, Cgen[a,i,b,j])
                end
            end
        end
    end
    return (Cgen=Cgenod, Cnat=Cnatod, Pigen=Pigen, Pinat=Pinat)
    
end


function pseudocount1(f1, pc, q::Int)
     return ((1-pc) .* f1 ) .+ (pc / q)
end

function pseudocount2(f2, pc, q::Int)
    return ((1-pc) .* f2 ) .+ (pc / q^2)
end


function reweight(msa,theta)
    L,M = size(msa)
    res = zeros(M)
    for m in 1:1000#M
        res[m] = 1 / sum(ham_dist(msa[:,m],msa)./L  .< theta)
    end
    return res
end

function conn_corr(msa, nat_msa, q::Int, L::Int; nat_weight::Float64 = 0.2, sim_weight::Float64 = 0.)

    f1_nat_ql, f2_nat_ql = compute_weighted_frequencies(nat_msa, 22, nat_weight); 
    f1_nat = reshape(f1_nat_ql,q,L); 
    f2_nat = reshape(f2_nat_ql,q,L,q,L);
    
    f1_ql, f2_ql = compute_weighted_frequencies(msa, 22, sim_weight); 
    f1 = reshape(f1_ql,q,L); 
    f2 = reshape(f2_ql,q,L,q,L);   

    c_ql = triu(f2_ql - f1_ql * f1_ql', 21); c_nat_ql = triu(f2_nat_ql - f1_nat_ql * f1_nat_ql', 21);

    c_nat = []; c = []; 
    for i in 1:L
        for j in i+1:L
            for a in 1:q
                for b in 1:q
                    push!(c_nat, f2_nat[a,i,b,j] - f1_nat[a,i] * f1_nat[b,j])
                    push!(c, f2[a,i,b,j] - f1[a,i] * f1[b,j])                    
                end
            end
        end
    end


    println(cor(c_nat[:], c[:]))
    println(cor(c_nat_ql[:], c_ql[:]))
    
end



function conn_corr(ff1::Array{T,1}, ff2::Array{T,1}, ff1_nat::Array{T,1}, ff2_nat::Array{T,2}) where {T}
 
    ql = length(ff1);
    q = 21; L = round(Int, ql/q);
    
    f1_nat = reshape(ff1_nat,q,L); 
    f2_nat = reshape(ff2_nat,q,L,q,L);
    
    f1 = reshape(ff1,q,L); 
    f2 = reshape(ff2,q,L,q,L);   

    c_nat = []; c = []; 
    for i in 1:L
        for j in i+1:L
            for a in 1:q
                for b in 1:q
                    push!(c_nat, f2_nat[a,i,b,j] - f1_nat[a,i] * f1_nat[b,j])
                    push!(c, f2[a,i,b,j] - f1[a,i] * f1[b,j])                    
                end
            end
        end
    end


    return c, c_nat
    
end


function conn_corr(ff1::Array{T,1}, ff2::Array{T,2}) where {T}
 
    ql = length(ff1);
    q = 21; L = round(Int, ql/q);
    
    cc = reshape(ff2 - ff1*ff1', q, L, q, L)
    c = []; 
    for i in 1:L
        for j in i+1:L
            for a in 1:q
                for b in 1:q
                    push!(c, cc[a,i,b,j])                    
                end
            end
        end
    end
    
    return c
end


function conn_corr_with_weight(msa, nat_msa, q, L, w)

    
    f1_nat_ql, f2_nat_ql = compute_weighted_frequencies(nat_msa, w, 22); 
    f1_nat = reshape(f1_nat_ql,q,L); 
    f2_nat = reshape(f2_nat_ql,q,L,q,L);
    
    f1_ql, f2_ql = compute_weighted_frequencies(msa, 22, 0.); 
    f1 = reshape(f1_ql,q,L); 
    f2 = reshape(f2_ql,q,L,q,L);   

    c_ql = triu(f2_ql - f1_ql * f1_ql', 21); c_nat_ql = triu(f2_nat_ql - f1_nat_ql * f1_nat_ql', 21);

    c_nat = []; c = []; 
    for i in 1:L
        for j in i+1:L
            for a in 1:q
                for b in 1:q
                    push!(c_nat, f2_nat[a,i,b,j] - f1_nat[a,i] * f1_nat[b,j])
                    push!(c, f2[a,i,b,j] - f1[a,i] * f1[b,j])                    
                end
            end
        end
    end


    println(cor(c_nat[:], c[:]))
    println(cor(c_nat_ql[:], c_ql[:]))
    
end


function conn_corr_with_pc(msa, nat_msa, q, L; pc = 10^-6, nat_weight = 0.2, sim_weight = 0.)

    f1_nat_ql, f2_nat_ql = compute_weighted_frequencies(nat_msa, 22, nat_weight); f1_nat = reshape(f1_nat_ql,q,L); f2_nat = reshape(f2_nat_ql,q,L,q,L);
    f1_nat_ql = pseudocount1(f1_nat_ql, pc, 21);
    f2_nat_ql = pseudocount2(f2_nat_ql, pc, 21);    
    f1_nat = pseudocount1(f1_nat, pc, 21);
    f2_nat = pseudocount2(f2_nat, pc, 21);
    
    f1_ql, f2_ql = compute_weighted_frequencies(msa, 22, sim_weight); f1 = reshape(f1_ql,q,L); f2 = reshape(f2_ql,q,L,q,L);
    f1_ql = pseudocount1(f1_ql, pc, 21);
    f2_ql = pseudocount2(f2_ql, pc, 21);
    f1 = pseudocount1(f1, pc, 21);
    f2 = pseudocount2(f2, pc, 21);
    
    c_ql = triu(f2_ql - f1_ql * f1_ql', 21); c_nat_ql = triu(f2_nat_ql - f1_nat_ql * f1_nat_ql', 21);

    c_nat = []; c = []; 
    for i in 1:L
        for j in 1:L
            if i !== j
                for a in 1:q
                    for b in 1:q
                        push!(c_nat, f2_nat[a,i,b,j] - f1_nat[a,i] * f1_nat[b,j])
                        push!(c, f2[a,i,b,j] - f1[a,i] * f1[b,j])
                    end
                end
            end
        end
    end


    println(cor(c_nat[:], c[:]))
    println(cor(c_nat_ql[:], c_ql[:]))
    
end









function check_equilibration(folder::String, nat_msa::Array{Int8,2}, step_msa::Array{Array{Int8,2},1}, steps::Array{Int,1}; nat_weight = 0.2, sim_weight = 0.)
    
    N_msa = length(steps)
    L, M = size(step_msa[end])
    
    if size(nat_msa,1) !== L
        error("Length of natural and artificial sequences is different \n
            remember that you need to have the msa in a (L, M) format of Int8")
    elseif N_msa !== length(step_msa)
        error("Length of step_msa and of steps is different")
    end
    
    f1_nat,f2_nat = compute_weighted_frequencies(nat_msa, 22, nat_weight); conn_nat = triu(f2_nat - f1_nat * f1_nat', 21);
    w = compute_weights(nat_msa, 22, nat_weight)[1];
    dist_nat = mean(ham_dist(step_msa[1][:,1], nat_msa), weights(w))
    err_dist_nat = std(ham_dist(step_msa[1][:,1], nat_msa), weights(w))
    cor1 = []; cor2 = []; corconn = []; dist = []; err_dist = [];
    a = zeros(21*L,21*L);
    for i in 1:N_msa
        f1,f2 = compute_weighted_frequencies(step_msa[i], 22, sim_weight); 
        conn = triu(f2 - f1 * f1',21);
        if i == N_msa
            a .= conn
        end
        push!(cor1, cor(f1[:], f1_nat[:]))
        push!(cor2, cor(f2[:], f2_nat[:]))
        push!(corconn, cor(conn[:], conn_nat[:]))
        dists = ham_dist(step_msa[1], step_msa[i])       
        push!(dist, mean(dists))
        push!(err_dist, var(dists))
    end 
    
    f1_end,f2_end = DCAUtils.compute_weighted_frequencies(step_msa[end], 22, sim_weight); 
    
    pc_nat, PC, expl_var = perform_pca(nat_msa, n_dim = 2);
    pc_sil = get_projection(PC, step_msa[end]);
            
    close("all")
    plt.plot(steps, cor1, label = "1-point"); plt.plot(steps, cor2, label = "2-point"); plt.plot(steps, corconn, label = "Conn cor"); plt.legend(); plt.xlabel("MCMC steps"); plt.ylabel("Pearson correlation"); plt.xscale("log");savefig(joinpath(folder,"freq_corr_evol.png"));
    
    
    close("all")
    plt.plot(f1_nat[:], f1_end[:], "o", label = "1-point"); plt.legend(); plt.xlabel("Natural frequencies"); plt.ylabel("Artificial frequencies"); savefig(joinpath(folder,"final_freq_natvssim.png"));
    
    
    close("all")
    plt.plot(steps, dist./L, label = "Simulated"); plt.plot(steps, [dist_nat for _ in 1:N_msa] ./L, label = "Natural"); plt.legend(); plt.xlabel("MCMC steps"); plt.ylabel("Mean hamming from wt (%)"); plt.xscale("log");savefig(joinpath(folder,"mean_hamming.png"))

    
    close("all")
    plt.plot(steps, (err_dist)./L, label = "Simulated"); plt.plot(steps, [err_dist_nat^2 for _ in 1:N_msa] ./L, label = "Natural"); plt.legend(); plt.xlabel("MCMC steps"); plt.ylabel("Variance hamming from wt (%)");plt.xscale("log"); savefig(joinpath(folder,"var_hamming.png"))
    
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
    

function check_equilibration(folder::String, nat_msa_path::String, step_msa::Array{Array{Int8,2},1}, steps::Array{Int,1})
    
    nat_msa = read_fasta_alignment(nat_msa_path, 0.9)
    check_equilibration(folder, nat_msa, step_msa, steps)
end


function check_equilibration_with_weight(folder::String, nat_msa::Array{Int8,2}, step_msa::Array{Array{Int8,2},1}, steps::Array{Int,1}, w; nat_weight = 0.2, sim_weight = 0.)
    
    N_msa = length(steps)
    L, M = size(step_msa[end])
    
    if size(nat_msa,1) !== L
        error("Length of natural and artificial sequences is different \n
            remember that you need to have the msa in a (L, M) format of Int8")
    elseif N_msa !== length(step_msa)
        error("Length of step_msa and of steps is different")
    end
    
    f1_nat,f2_nat = compute_weighted_frequencies(nat_msa, w, 22); conn_nat = triu(f2_nat - f1_nat * f1_nat', 21);
    #w = compute_weights(nat_msa, 22, nat_weight)[1];
    dist_nat = mean(ham_dist(step_msa[1][:,1], nat_msa), weights(w))
    err_dist_nat = std(ham_dist(step_msa[1][:,1], nat_msa), weights(w))
    cor1 = []; cor2 = []; corconn = []; dist = []; err_dist = [];
    a = zeros(21*L,21*L);
    for i in 1:N_msa
        f1,f2 = compute_weighted_frequencies(step_msa[i], 22, sim_weight); 
        conn = triu(f2 - f1 * f1',21);
        if i == N_msa
            a .= conn
        end
        push!(cor1, cor(f1[:], f1_nat[:]))
        push!(cor2, cor(f2[:], f2_nat[:]))
        push!(corconn, cor(conn[:], conn_nat[:]))
        dists = ham_dist(step_msa[1], step_msa[i])       
        push!(dist, mean(dists))
        push!(err_dist, var(dists))
    end 
    
    f1_end,f2_end = DCAUtils.compute_weighted_frequencies(step_msa[end], 22, sim_weight); 
    
    pc_nat, PC, expl_var = perform_pca(nat_msa, n_dim = 2);
    pc_sil = get_projection(PC, step_msa[end]);
            
    close("all")
    plt.plot(steps, cor1, label = "1-point"); plt.plot(steps, cor2, label = "2-point"); plt.plot(steps, corconn, label = "Conn cor"); plt.legend(); plt.xlabel("MCMC steps"); plt.ylabel("Pearson correlation"); plt.xscale("log");savefig(joinpath(folder,"freq_corr_evol.png"));
    
    
    close("all")
    plt.plot(f1_nat[:], f1_end[:], "o", label = "1-point"); plt.legend(); plt.xlabel("Natural frequencies"); plt.ylabel("Artificial frequencies"); savefig(joinpath(folder,"final_freq_natvssim.png"));
    
    
    close("all")
    plt.plot(steps, dist./L, label = "Simulated"); plt.plot(steps, [dist_nat for _ in 1:N_msa] ./L, label = "Natural"); plt.legend(); plt.xlabel("MCMC steps"); plt.ylabel("Mean hamming from wt (%)"); plt.xscale("log");savefig(joinpath(folder,"mean_hamming.png"))

    
    close("all")
    plt.plot(steps, (err_dist)./L, label = "Simulated"); plt.plot(steps, [err_dist_nat^2 for _ in 1:N_msa] ./L, label = "Natural"); plt.legend(); plt.xlabel("MCMC steps"); plt.ylabel("Variance hamming from wt (%)");plt.xscale("log"); savefig(joinpath(folder,"var_hamming.png"))
    
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
    

function check_equilibration_with_weight(folder::String, nat_msa_path::String, step_msa::Array{Array{Int8,2},1}, steps::Array{Int,1}, w)
    
    nat_msa = read_fasta_alignment(nat_msa_path, 0.9)
    check_equilibration_with_weight(folder, nat_msa, step_msa, steps, w)
end