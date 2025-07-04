function analyze_msa(msa, target_msa, wt, h, J, round; n_seq = 100)
    L = length(wt)
    target_entr = get_entropy(target_msa, q = 20)
    entr = get_entropy(msa,q=20);
    corr = cor(entr, target_entr); 
    en = energy(msa,h,J)
    ham = ham_dist(wt,msa) ./L
    pair_ham = pairwise_ham_dist(msa, n_seq = n_seq, all = true) ./L
    
    return (entr = entr, en = en, ham = ham, pair_ham = pair_ham, corr = corr, round = round) 
end

function analyze_target_msa(msa, wt, h, J, round; n_seq = 100)
    L = length(wt)
    entr = get_entropy(msa,q=20);
    en = energy(msa,h,J)
    ham = ham_dist(wt,msa) ./L
    pair_ham = pairwise_ham_dist(msa, n_seq = n_seq, all = true) ./L
    
    return (entr = entr, en = en, ham = ham, pair_ham = pair_ham) 
end

    
    
    
    
    