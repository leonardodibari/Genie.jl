#export fasta2matrix, map_gaps, read_par_BM, read_par_BM_LR

id(i, a, q) = (i .- 1).*q .+ a

function read_graph_new(filepath)
    file = open(filepath, "r"); content = read(file, String)[1:end-1]; close(file)
        content = split(content, "\n")
            alphabet, idx = "", 0
            for line ∈ content
                if string(line[1]) == "h"
                    _, i, a, value = split(line, " ")
                    if alphabet == "" || !occursin(a, string(alphabet))
                        alphabet = alphabet * a
                    end
                    idx += 1
                end
            end
            Nq = length(alphabet)
            Nv = div(idx, Nq)
            J, vbias = zeros(Float32, Nv*Nq, Nv*Nq), zeros(Float32, Nq, Nv)
            for line ∈ content
                if string(line[1]) == "J"
                    _, i, j, a, b, value = split(line, " ")
                    i, j = parse(Int64, i) + 1, parse(Int64, j) + 1
                    id1, id2 = id(i, findfirst(a, string(alphabet)), Nq)[1], id(j, findfirst(b, string(alphabet)),Nq)[1]
                    J[id1, id2], J[id2, id1] = parse(Float32, value), parse(Float32, value)
                elseif string(line[1]) == "h"
                    _, i, a, value = split(line, " ")
                    i = parse(Int64, i) + 1
                    vbias[findfirst(a, string(alphabet))[1], i] = parse(Float32, value)
                end
            end
        return J, vbias, alphabet
    end


function read_par_BM_standardorder(path::AbstractString, q::Integer = 21)
    params = readdlm(path,' ', use_mmap = true)[:, 2:6]
    l_file = size(params, 1)
    N = Integer(((q - 2) + sqrt( (q-2)^2 + 8*l_file))/(2*q))
    J = Array{Float64}(undef, q, q, N, N)
    h = Array{Float64}(undef, q, N)
    n_J = Int(q*q*N*(N-1)/2)
    n_h = q*N
    for k in 1:n_J
        i, j, a, b, par_j = params[k, :]
        i += 1
        j += 1
        a += 1
        b += 1
        J[a, b, i, j] = par_j
    end
    for l in (n_J + 1): n_h + n_J
        i, a, par_h = params[l, :]
        i += 1
        a += 1
        h[a, i] = par_h
    end
    return h, J
end

"""
        function map_gaps(matrix)
        
        Returns a matrix with 21's instead of 0's,
        to change between different mappings.        

"""

function map_gaps(matrix)
        for i in 1:size(matrix)[1]
		for j in 1:size(matrix)[2]
			if (matrix[i, j] == 0) 
				matrix[i, j] = 21
			end
		end
	end
	return matrix
end



export read_par_BM

"""
	read_par_BM(path::AbstractString, q, N)

	Reads the parameters of a Potts model in format
	
	J i j a b
	...
	h i a 

	and returns them in tensor format.
	J[a, b, i, j] is such that  1 == "A", ... "21" == "-",
	and the same for  h[a, i].

"""

function read_par_BM(path::AbstractString, q, N)
	data = readdlm(path,' ', use_mmap = true)[:, 2:6]
	J = Array{Float64}(undef, q, q, N, N)
	h = Array{Float64}(undef, q, N)
	n_J = Int(q*q*N*(N-1)/2)
	n_h = q*N

	for k in 1:n_J
		i, j, a, b, par_j = data[k, :]
		i += 1
		j += 1
		a == 0 && (a = 21)
		b == 0 && (b = 21)
		J[a, b, i, j] = par_j
	end

	for l in (n_J + 1): n_h + n_J
		i, a, par_h = data[l, :]
		i += 1
		a == 0 && (a = 21)
		h[a, i] = par_h
	end

	return h, J
end






"""
    vec2string(v)

    Takes as input a vector of integers (representing an amino acids)
    and returns a list of characters, the corresponding amino acids. 
    In this case with the convention "1 2 .. 21" == "A C .. -"
"""


function vec2string(v)
    s = ""
    for i in v
        s = s*num2letter(i)
    end
    return s
end


"""
    manage_fields(h)

    Takes in input a 2-dimensional tensor, a matrix,
    containing the values of the fiels for every amino acid a
    and position i. 
    Returns the same matrix, by removing the biggerst field value
    for every position i, so that the biggest field results 0.
    Help in computing exponentials without going to overflow.

"""


function manage_fields(h_old)
   h_new = deepcopy(h_old)
   q, N = size(h_old)
   for i in 1:N
        hmax = maximum([h_old[a, i] for a in 1:21])
        h_new[:, i] .-= hmax
   end
   return h_new
end



"""
    simmetrize_couplings(J)

    Takes in input a 4-dimensional tensor,
    containing the values of the couplings for every couple of 
    amino acids a, b and positions i, j. 
    Returns the same tensor, by simmetrizing it, 
    since J[a, b, i, j] = 0 if i < j, or viceversa, depends.

"""

function simmetrize_J(J_old)
   J_s = deepcopy(J_old)
   q,q, N, N = size(J_old)
   for i in 1:N
        for j in i+1:N
            for a in 1:q
                for b in 1:q
                      J_s[b, a, j, i]  = J_old[a, b, i, j] 
                end
            end
        end
    end
    return J_s
end


"""
    read_DNA(file_path::AbrastractString)

    Takes path of file with dna sequences and extracts a msa of dna sequences
"""
function read_DNA(file_path::AbstractString)
    # Initialize an empty vector to store sequences
    sequences = Vector{String}()
    
    # Read the file
    open(file_path) do file
        # Initialize a variable to store the current sequence
        current_sequence = ""
        
        for line in eachline(file)
            # Check if the line is empty or a header line
            if isempty(line) || startswith(line, '>')
                # If the current sequence is not empty, store it
                if !isempty(current_sequence)
                    push!(sequences, current_sequence)
                    current_sequence = ""
                end
                continue
            end
            
            # Add the line to the current sequence
            current_sequence *= strip(line)
        end
        
        # Store the last sequence
        if !isempty(current_sequence)
            push!(sequences, current_sequence)
        end
    end
    
    # Determine the maximum length of sequences
    max_sequence_length = maximum(length.(sequences))
    
    # Initialize the matrix with the correct dimensions
    sequences_matrix = fill("", length(sequences), max_sequence_length ÷ 3)
    
    # Split each sequence into chunks of three characters and store in the matrix
    for (i, sequence) in enumerate(sequences)
        for (j, chunk) in enumerate([sequence[k:k+2] for k in 1:3:length(sequence)-2])
            sequences_matrix[i, j] = chunk
        end
    end
    
    return sequences_matrix
end


function build_dna_seq_steps(step_matrices::Vector{Matrix{String}})
    seq_matrices = Matrix{String}[]
    for i in 1:size(step_matrices[1], 1)
        seq_matrix = Matrix{String}(undef, length(step_matrices), size(step_matrices[1], 2))
        for j in 1:length(step_matrices)
            seq_matrix[j, :] = step_matrices[j][i, :]
        end
        push!(seq_matrices, seq_matrix)
    end
    return seq_matrices
end

function matrixtofasta(path, MSA)
    n_seq = size(MSA, 1)
    FastaWriter(path, "w") do file
        for i in 1:n_seq	
            seq = MSA[i, :]
            writeentry(file, "$i", vec2string(seq) )
        end
    end
end

function leafstofasta(path, tree)
    n_seq = length(tree.lleaves)
    FastaWriter(path, "w") do file
        for a in keys(tree.lleaves) 
            writeentry(file, "$(a)", vec2string(tree["$(a)"].data.seq[:]))
        end
    end
end


function nodestofasta(path, tree)
    n_seq = length(tree.lleaves)
    FastaWriter(path, "w") do file
        for a in keys(tree.lnodes) 
            writeentry(file, "$(a)", vec2string(tree["$(a)"].data.seq[:]))
        end
    end
end
