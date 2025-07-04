# translate ordered list of residues to paths of sequences
function path2seq(seq_start::Array{T, 1}, seq_end::Array{T, 1}, path_muts::Array{Int16, 1}) where {T}
    
    M = length(path_muts)
    L = length(seq_start)

    path_seq = T.(zeros(L, M+1))
    path_seq[:,1] = seq_start
    
    for i in 1:M
        ss = copy(path_seq[:,i])
        res = path_muts[i]
        ss[res] = seq_end[res]
        path_seq[:, i+1] = ss
    end
    
    return path_seq   
end


function delta_energy(h::Array{Float64,2},
                              J::Array{Float64,4},
                              S::Array{<:Integer,1}, 
                              ref::Array{<:Integer, 1})
    
    q, N = size(h)
    E = 0.0
    index_v = collect(1:N)
    bool_common = Array{Bool}(undef, N)
    ll = 0
    @fastmath for val in S
        ll+=1
        @inbounds bool_common[ll] = (val == ref[ll])
    end
    common = index_v[bool_common]
    non_common = index_v[(!).(bool_common)]
    
    @fastmath for i in non_common
         @inbounds  E -= (h[S[i],i] - h[ref[i],i])
         for j = 1:N
            if j > i
              @inbounds  E -= (J[S[j],j,S[i],i] - J[ref[j],j,ref[i],i] )
            end
        end
    end
    
    @fastmath for i in common
         for j in non_common
            if j > i
              @inbounds  E -= (J[S[j],j,S[i],i] - J[ref[j],j,ref[i],i] )
            end
        end
    end
    
    return E
end

# save to file non aligned sequences adding VIM2 terminal
function dealign_and_save_VIM2_NDM1_intermediates(path_in)
    
    VIM2_Nterm = "mfkllskllvyltasimaiasplafsvdssgeyptvseipv"
    VIM2_Cterm = "tnrsvve"

    MSA = fasta2matrix(path_in)
    M, L = size(MSA)
    desc = readlines(path_in)[1:4:end]
    path_out = path_in * "_full_sequence"
    
    open(path_out, "w") do file
        for m in 1:M
            seq = MSA[m, :]
            seq_string = join(num2letter.(seq))
            #seq_nog = join([seq_string[i] for i in 1:L if seq_string[i] != "-" ])
            seq_chimera = VIM2_Nterm * seq_string * VIM2_Cterm
            println(file, ">$(m) $(desc[m][4:end])")
            println(file, seq_chimera)
        end
    end
end

function meanE_path(seq_start::Array{T, 1}, seq_end::Array{T, 1},       
        h::Array{Float64, 2}, J::Array{Float64, 4};  β::Real = 1, steps::Int64 = 100) where {T}
    
    # set of non common positions (including gaps)
    not_common = [i for (i, amino) in enumerate(seq_end) if seq_start[i] != amino]
    
    # length of non common position
    L_not_common = length(not_common)

    # initialize path
    old_path_muts = Int16.(sample(not_common, L_not_common, replace = false))
    
    # index for iterations
    k = 0
    
    score = []
    
    while k < steps
        
        # add 1 to iteration index
        k += 1
        
        # new path
        new_path_muts = copy(old_path_muts)   
        
        # mutations to swap
        pos_res1, pos_res2 = sample(1:L_not_common, 2, replace = false)
        
        # introduce mutations in new path
        new_path_muts[pos_res1] = old_path_muts[pos_res2]
        new_path_muts[pos_res2] = old_path_muts[pos_res1]
        
        # build new trajectory
        new_path_seq = path2seq(seq_start, seq_end, new_path_muts)
        
        # build trajectory from the mutation's order
        old_path_seq = path2seq(seq_start, seq_end, old_path_muts)

        # compute energy of new trajectory 
        delta_en = mean([delta_energy(h, J, new_path_seq[:, i], old_path_seq[:, i]) for i in 1:(L_not_common+1) ])

        # metropolis acceptance rule
        if min(1, exp( -β* (delta_en) )) > rand()
            old_path_muts = new_path_muts
        end   
        push!(score, mean(energy(new_path_seq, h, J)))
    end
    
    return old_path_muts, score
end


function meanE_path_ij(pos1::Int, pos2::Int, seq_start::Array{T, 1}, seq_end::Array{T, 1},       
        h::Array{Float64, 2}, J::Array{Float64, 4}; q = 21, reg = 0.01, β::Real = 1, steps::Int64 = 100) where {T}
    
    
    score = []
    L = length(seq_start)
    
    # set of non common positions (including gaps)
    not_common = [i for (i, amino) in enumerate(seq_end) if seq_start[i] != amino]
    
    filter!(x -> x !== pos1,  not_common)
    filter!(x -> x !== pos2,  not_common)
    filter!(x -> seq_start[x] !== 21,  not_common)
    filter!(x -> seq_end[x] !== 21,  not_common)
    
    # length of non common position
    L_not_common = length(not_common)

    # initialize path
    old_path_muts = Int16.(sample(not_common, L_not_common, replace = false))
    
    # index for iterations
    k = 0
    
    while k < steps
        
        # add 1 to iteration index
        k += 1
        
        # new path
        new_path_muts = copy(old_path_muts)   
        
        # mutations to swap
        pos_res1, pos_res2 = sample(1:L_not_common, 2, replace = false)
        
        # introduce mutations in new path
        new_path_muts[pos_res1] = old_path_muts[pos_res2]
        new_path_muts[pos_res2] = old_path_muts[pos_res1]
        
        # build new trajectory
        new_path_seq = path2seq(seq_start, seq_end, new_path_muts)
        
        # build trajectory from the mutation's order
        old_path_seq = path2seq(seq_start, seq_end, old_path_muts)

        # compute energy of new trajectory 
        delta_en = mean([delta_energy(h, J, new_path_seq[:, i], old_path_seq[:, i]) for i in 1:(L_not_common+1) ])
        
        p_ij_old = [dist_proba_2_point(pos1, pos2, old_path_seq[:,n], h, J, L, q, 1.)[seq_end[pos1], seq_end[pos2]] for n in 1:(L_not_common+1)]
        
        p_ij_new = [dist_proba_2_point(pos1, pos2, new_path_seq[:,n], h, J, L, q, 1.)[seq_end[pos1], seq_end[pos2]] for n in 1:(L_not_common+1)]
        
        
        delta_p_ij = mean(p_ij_new - p_ij_old)

        # metropolis acceptance rule
        if min(1, exp( -β* (-delta_p_ij + reg*delta_en) )) > rand()
            old_path_muts = new_path_muts
        end   
        push!(score, mean(p_ij_new))
    end
    
    return old_path_muts, score
end


function sample_n_paths(seq1, seq2, h::Array{Float64, 2}, J::Array{Float64, 4}, npath::Int; n_steps = 5000, β = 5)

    # non common residue between VIM2 and NDM1
    not_common_residues = [i for (i, amino) in enumerate(seq1) if seq2[i] != amino]
    L = length(not_common_residues)

    # initialize 3 paths
    mutations_rand = Int16.(sample(not_common_residues, L, replace = false))
    paths_meanE_vec = [path2seq(seq1, seq2, mutations_rand) for i in 1:npath]

    # initialize energy of the path
    e_meanE_vec = [zeros(L + 1) for i in 1:npath]
    scores = []
    # generate 2 paths
    for i in 1:npath
        @time final_muts, score = meanE_path(seq1, seq2, h, J, β = β, steps = n_steps) 
        paths_meanE_vec[i] = path2seq(seq1, seq2, final_muts)
        e_meanE_vec[i] = energy(paths_meanE_vec[i], h, J)
        push!(scores, score)
    end
    return (seqs = paths_meanE_vec, ens = e_meanE_vec, scores = scores)
end



function sample_n_paths_ij(pos1::Int, pos2::Int, sseq1, seq2, h::Array{Float64, 2}, J::Array{Float64, 4}, npath::Int; reg = 0.01, n_steps = 5000, β = 5, q = 21)

    
    L_tot = length(sseq1)
    
    seq1 = deepcopy(sseq1);
    seq1[pos1] = seq2[pos1]; seq1[pos2] = seq2[pos2];
    not_common_residues = findall(seq1 .!= seq2) 
    filter!(x -> x !== pos1,  not_common_residues)
    filter!(x -> x !== pos2,  not_common_residues)
    filter!(x -> seq1[x] !== 21,  not_common_residues)
    filter!(x -> seq2[x] !== 21,  not_common_residues)
    
    L = length(not_common_residues)

    # initialize paths
    mutations_rand = Int16.(sample(not_common_residues, L, replace = false))
    paths_meanE_vec = [path2seq(seq1, seq2, mutations_rand) for i in 1:npath]

    # initialize energy of the path
    e_meanE_vec = [zeros(L + 1) for i in 1:npath]
    ps_ij = [zeros(L + 1) for i in 1:npath]
    scores = [zeros(L + 1) for i in 1:npath]
    muts = [zeros(L + 1) for i in 1:npath]

    # generate 2 paths
    @tasks for i in 1:npath
        final_muts, score = meanE_path_ij(pos1, pos2, seq1, seq2, h, J, β = β, reg = reg,  steps = n_steps) 
        paths_meanE_vec[i] = path2seq(seq1, seq2, final_muts)
        e_meanE_vec[i] = energy(paths_meanE_vec[i], h, J)
        ps_ij[i] = [dist_proba_2_point(pos1, pos2, paths_meanE_vec[i][:,n], h, J, L_tot, q, 1.)[seq2[pos1], seq2[pos2]] for n in 1:size(paths_meanE_vec[i],2)]
        scores[i] = score
        muts[i] = final_muts
    end
    
    return (seqs = paths_meanE_vec, ens = e_meanE_vec, ps_ij = ps_ij, scores = scores, muts = muts)
end





function pairwise_hd(MSA)
    M, L = size(MSA)
    HD = zeros( Int64.(M * (M-1) / 2))
    k = 0
    for i in 1:(M-1)
        for j in (i+1):M
            k += 1
            HD[k] = ham_dist(MSA[:,i], MSA[:, j])
        end
    end
    return HD
end

