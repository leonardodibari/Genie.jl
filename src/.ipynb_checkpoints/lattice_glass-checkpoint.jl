function idx_to_ijk(idx::Int, L::Int)
    
    k = div(idx - 1, m * n) + 1
    j = div(idx - 1, m) % n + 1
    i = (idx - 1) % m + 1
    
    return [i, j, k]
end

function swap_occupation!(lattice::Array{Int,3}, pos1::Array{Int,1}, pos2::Array{Int,1})
    part_type = lattice[pos1[1], pos1[2], pos1[3]]
    lattice[pos1[1], pos1[2], pos1[3]] = lattice[pos2[1], pos2[2], pos2[3]] 
    lattice[pos2[1], pos2[2], pos2[3]] = part_type ;
end
 


function initialize_lattice(frac_A::Float64, frac_free::Float64, L::Int)
    
    lattice = zeros(Int, L, L, L)

    total_sites = L^3

    # Calculate number of free sites and particle counts
    num_free = Int(frac_free * total_sites)
    num_A = Int(frac_A * (total_sites - num_free))
    num_B = total_sites - num_free - num_A

    # Create a list of all lattice positions
    positions = [(x, y, z) for x in 1:L, y in 1:L, z in 1:L]
    shuffle!(positions)  # Shuffle the positions randomly

    # Place particles of type A
    for pos in positions[1:num_A]
        lattice[pos...] = 1
    end

    # Place particles of type B
    for pos in positions[num_A+1:num_A+num_B]
        lattice[pos...] = 2
    end

    return lattice
end


function distance!(arr1::Array{Int,1}, arr2::Array{Int,1})
    return sqrt(sum(abs2, arr1 .- arr2))
end

function swap_options(pos::Array{Int,1}, lattice::Array{Int,3}, L::Int)
    x = pos[1]; y = pos[2]; z = pos[3];
    a = lattice[x,y,z]
    list = []
    if lattice[mod1(x + 1, L), y, z] != a push!(list, [mod1(x + 1, L), y, z]) end  # Right
    if lattice[mod1(x - 1, L), y, z] != a push!(list, [mod1(x - 1, L), y, z]) end  # Left
    if lattice[x, mod1(y + 1, L), z] != a push!(list, [x, mod1(y + 1, L), z]) end  # Up        
    if lattice[x, mod1(y - 1, L), z] != a push!(list, [x, mod1(y - 1, L), z]) end  # Down            
    if lattice[x, y, mod1(z + 1, L)] != a push!(list, [x, y, mod1(z + 1, L)]) end  # Forward
    if lattice[x, y, mod1(z - 1, L)] != a push!(list, [x, y, mod1(z - 1, L)]) end  # Backward
    
    return list 
end

function neighbours(pos::Array{Int,1}, lattice::Array{Int,3}, L::Int)
    x = pos[1]; y = pos[2]; z = pos[3];
    list = []
    if lattice[mod1(x + 1, L), y, z] != 0 push!(list, [mod1(x + 1, L), y, z]) end  # Right
    if lattice[mod1(x - 1, L), y, z] != 0 push!(list, [mod1(x - 1, L), y, z]) end  # Left
    if lattice[x, mod1(y + 1, L), z] != 0 push!(list, [x, mod1(y + 1, L), z]) end  # Up        
    if lattice[x, mod1(y - 1, L), z] != 0 push!(list, [x, mod1(y - 1, L), z]) end  # Down            
    if lattice[x, y, mod1(z + 1, L)] != 0 push!(list, [x, y, mod1(z + 1, L)]) end  # Forward
    if lattice[x, y, mod1(z - 1, L)] != 0 push!(list, [x, y, mod1(z - 1, L)]) end  # Backward
    
    return list 
end
    
function count_neighbours(pos::Array{Int,1}, lattice::Array{Int,3}, L::Int)    
    list = neighbours(pos, lattice, L)
    return length(list)
end
    
function node_energy(pos::Array{Int,1}, lattice::Array{Int,3}, pref::Array{Int,1}, L::Int)
    return (count_neighbours(pos, lattice, L) - pref[lattice[pos[1],pos[2],pos[3]]])^2
end
    

function select_pos(rng::Xoshiro, lattice::Array{Int,3}, L::Int)
    pos1 = rand(rng, 1:L,3)
    list = swap_options(pos1, lattice, L)
    if length(list) == 0
        return select_pos(rng::Xoshiro, lattice::Array{Int,3}, L::Int)
    else
        pos2 = list[rand(rng, 1:length(list))]
        return pos1,pos2
    end 
end
    

function lattice_delta_energy(pos1::Array{Int,1}, pos2::Array{Int,1}, lattice::Array{Int,3}, lattice_f::Array{Int,3}, pref::Array{Int,1}, L::Int)
    
    if lattice[pos1[1], pos1[2], pos1[3]] * lattice[pos2[1], pos2[2], pos2[3]] > 0
        n1 = count_neighbours(pos1, lattice, L)
        n2 = count_neighbours(pos2, lattice, L)
        l1 = pref[lattice[pos1[1], pos1[2], pos1[3]]]
        l2 = pref[lattice[pos2[1], pos2[2], pos2[3]]]
        return 2*(l1-l2)*(n1-n2)
    end
    
    en_i = 0; en_f = 0; 
    
    #lattice2 .= lattice
    lattice_f[pos1[1], pos1[2], pos1[3]] = lattice[pos2[1], pos2[2], pos2[3]] 
    lattice_f[pos2[1], pos2[2], pos2[3]] = lattice[pos1[1], pos1[2], pos1[3]]
    
    ne1_i = neighbours(pos1, lattice, L); 
    ne1_f =  neighbours(pos1, lattice_f, L);
    ne2_i = neighbours(pos2, lattice, L); 
    ne2_f =  neighbours(pos2, lattice_f, L);
    
    for p in ne1_i
        en_i += node_energy(p, lattice, pref, L)
    end 
    for p in ne1_f
        en_f += node_energy(p, lattice_f, pref, L)
    end
    
    for p in ne2_i
        en_i += node_energy(p, lattice, pref, L)
    end
    for p in ne2_f
        en_f += node_energy(p, lattice_f, pref, L)
    end
        
    return en_f - en_i   
end

function lattice_energy(lattice::Array{Int64, 3}, pref::Array{Int,1}, L::Int)
    en = 0
    
    # Iterate through each site in the lattice
    for x in 1:L, y in 1:L, z in 1:L
        # Skip empty sites (i.e., value 0)
        if lattice[x, y, z] !== 0          
            en += node_energy([x,y,z], lattice, pref, L)
        end
    end

    return en
end



function glass_metropolis!(pos1::Array{Int,1}, pos2::Array{Int,1}, lattice::Array{Int64, 3}, lattice2::Array{Int,3}, pref::Array{Int,1}, rng::Xoshiro, T::Float64, L::Int)
    
    if rand(rng) < exp(-lattice_delta_energy(pos1, pos2, lattice, lattice2, pref, L)/T)
        swap_occupation!(lattice, pos1, pos2)
    else
        swap_occupation!(lattice2, pos1, pos2)
    end
end
        

function run_lattice_dynamics(start_latt; 
        N_steps::Int = 100, 
        L::Int = 20,
        phi1::Float64 = 0.4,
        phi_free::Float64 = 0.25,
        pref::Array{Int,1} = [3,5],
        T::Float64 = 1.0,
        rand_init = true,
        N_points::Union{Int, Nothing} = nothing, 
        each_step::Union{Int, Nothing} = nothing, 
        verbose = true)
    
    N_chains = length(start_latt)
    L = size(start_latt[1],1)
    
    if N_points !== nothing && each_step !== nothing
        error("You cannot specify both N_points and each_step: 
            \n-if you specify N_points, you get some msa (N_points) in logarithmic scale along the trajectory
            \n-if you specify each_step, you get a msa every each_step along the trajectory
            \n-if you don't specify anything then you get only the msa at the final step \n")
    end
    
    rng = random_gens(N_chains)
    count = 0

    if rand_init == true
        println("Random Initialization")
        chains = [initialize_lattice(phi1, phi_free, L) for n in 1:N_chains]
        chains2 = deepcopy(chains)
    else
        chains = start_latt 
        chains2 = deepcopy(chains)
    end
     
    if N_points !== nothing 
        if N_points > N_steps
            error("N_points must be smaller than N_steps")
        end
        steps = unique([trunc(Int,10^y) for y in range(log10(1), log10(N_steps), length=N_points)])
        step_latt = zeros(Int, L^3, length(steps), N_chains)
    end
    
    if each_step !== nothing 
        if each_step > N_steps
            error("each_step must be smaller than N_steps")
        end
        steps = [i for i in 1:each_step:N_steps]
        step_latt = zeros(Int, L^3, length(steps), N_chains)
    end
    
     
    @inbounds for t in 1:N_steps
        
        if ((N_points !== nothing) || (each_step !== nothing)) && (t in steps)
            count += 1
            if verbose == true
                println(t)
            end
            @tasks for n in 1:N_chains
                for i in 1:L^3
                    step_latt[i,count,n] = chains[n][i]
                end
            end
        end
        
        @tasks for n in 1:N_chains
            pos1, pos2 = select_pos(rng[n], chains[n], L)
            glass_metropolis!(pos1, pos2, chains[n], chains2[n], pref, rng[n], T, L)            
        end
    end  
    
    if (N_points !== nothing) || (each_step !== nothing) 
        return (step_latt = step_latt, steps = steps, T = T)
    else
        lattice = zeros(Int, L^3, N_chains)
        @tasks for n in 1:N_chains
            for i in 1:L^3
                lattice[i,n] = chains[n][i]
            end
        end
        return (lattice = lattice, T = T)
    end 
end

function copy_lattice!(dest::Array{Int,3}, source::Array{Int,3}, L::Int)
    for i in 1:L
        for j in 1:L
            for k in 1:L
                dest[i,j,k] = source[i,j,k]
            end
        end
    end
end
    

function lattice_ham(latt1::Array{Int,1}, latt2::Array{Int,1}, L::Int)
    return sum(latt1 .!= latt2)
end

function lattice_ham_on_chain(step_latt::Array{Int,2})
    LL, N_points = size(step_latt)
    L = round(Int,LL^(1/3))
    res = zeros(N_points)
    for n in 1:N_points
        res[n] = lattice_ham(step_latt[:,1], step_latt[:,n],L)
    end
    return res
end


function all_ham_dist(step_latt::Array{Int,3})
    N_chains = size(step_latt,3)
    return Float64.(hcat([Genie.lattice_ham_on_chain(step_latt[:,:,i]) for i in 1:N_chains]...))
end

    
    
    
   

        
    
    
    