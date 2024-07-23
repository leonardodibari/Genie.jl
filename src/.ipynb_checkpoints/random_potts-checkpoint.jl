function run_potts(start_msa, h::Array{T,2}, J::Array{T,4}; 
        N_steps::Int = 100, 
        temp = 1.0,
        N_points::Union{Int, Nothing} = nothing, 
        each_step::Union{Int, Nothing} = nothing, 
        rand_init = false, 
        q = 21,
        verbose = false) where {T}
    
    
    L, N_chains = size(start_msa)
    if (size(J,1) !== size(J,3)) || (size(J,2) !== size(J,4))
        error("Size of J should be (q,L,q,L)")
    elseif (size(J,2) !== L) || (size(h,2) !== L)
        error("Length of sequences different from length of parameters")
    elseif N_points !== nothing && each_step !== nothing
        error("You cannot specify both N_points and each_step: 
            \n-if you specify N_points, you get some msa (N_points) in logarithmic scale along the trajectory
            \n-if you specify each_step, you get a msa every each_step along the trajectory
            \n-if you don't specify anything then you get only the msa at the final step \n")
    end
    
    rng = random_gens(N_chains)
    count = 0
    temp = T(temp)
    
    
    if rand_init == true
        println("Random Initialization")
        chains = [AminoChain(Int8.(rand(1:q, L)), q, rng[n]) for n in 1:N_chains]
    else
        chains = [AminoChain(start_msa[:,n], q, rng[n]) for n in 1:N_chains]
    end
     
    if N_points !== nothing 
        if N_points > N_steps
            error("N_points must be smaller than N_steps")
        end
        steps = unique([trunc(Int,10^y) for y in range(log10(1), log10(N_steps), length=N_points)])
        step_msa = [zeros(Int8, (L, N_chains)) for i in 1:length(steps)]
    end
    
    if each_step !== nothing 
        if each_step > N_steps
            error("each_step must be smaller than N_steps")
        end
        steps = [i for i in 1:each_step:N_steps]
        step_msa = [zeros(Int8, (L, N_chains)) for i in 1:length(steps)]
    end
    
     
    @inbounds for t in 1:N_steps
        
        if ((N_points !== nothing) || (each_step !== nothing)) && (t in steps)
            count += 1
            if verbose == true
                println(t)
            end
            @tasks for n in 1:N_chains
                for i in 1:L
                    step_msa[count][i,n] = chains[n].seq[i]
                end
            end
        end
        
        run_simple_metropolis!(chains, h, J, N_chains, temp, L, q)
            
    end  
    
    if (N_points !== nothing) || (each_step !== nothing) 
        return (step_msa = step_msa, steps = steps, temp = temp)
    else
        msa = Int8.(zeros(L, N_chains))
        @tasks for n in 1:N_chains
            for i in 1:L
                msa[i,n] = Int8.(chains[n].seq[i])
            end
        end
        return (msa = msa, temp = temp)
    end 
end
    

function run_simple_metropolis!(chains, 
        h::Array{T,2}, 
        J::Array{T,4},
        N_chains::Int,
        temp::T,
        L::Int, 
        q::Int) where {T}
    
    @tasks for n in 1:N_chains
        
        seq_site = rand(chains[n].generator, 1:L)
        new_amino = rand(chains[n].generator, 1:q) 
        
        dE = single_mut_dE(chains[n].seq, h, J, new_amino, seq_site, L)
        if rand(chains[n].generator) < exp(-dE/temp)
            chains[n].seq[seq_site] = new_amino
        end
    end
    
end


function Potts_fixed_mean_conn(L::Int, C::Int, var_J::T, var_h::T, q::Int; mean_h::T = zero(T), mean_J::T = zero(T)) where {T}
    
    p = C/(L-1)
    
    ### Initialize seed RNG
    RNGseed = 10
    Random.seed!(RNGseed)


    ### Initialize fields & couplings
    FieldDistribution = Laplace(T(mean_h),var_h) 
    h=zeros(Float64,(q,L))
    
    for i in 1:L
        for a in 1:q
            h[a,i] = rand(FieldDistribution)
        end
    end

    J = zeros(Float64, (q,L,q,L))
    for i in 1:L
        for j in i+1:L
            if rand() < p
                for a in 1:q
                    for b in 1:q
                        Jvalue = rand(Normal(T(mean_J),var_J))
                        #Jvalue = 1
                        J[a,i,b,j] = Jvalue
                        J[b,j,a,i] = Jvalue
                        #=if a == b
                            J[a,i,b,j] = Jvalue
                            J[b,j,a,i] = Jvalue
                    end=#
                    end
                end
            end
        end
    end

    for i in 1:L
        for j in 1:L
            for a in 1:q
                for b in 1:q
                    @assert J[a,i,b,j] == J[b,j,a,i]
                end
            end
        end
    end
    
    return h,J
end

