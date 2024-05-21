module Assignment07

export normalizeDNA,       
        composition,
        gc_content,
        complement,
        reverse_complement,
        parse_fasta 

# # uncomment the following line if you intend to use BioSequences types
# using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end


"""
    composition(sequence)

Counts the number of each type of base
in a DNA sequence and returns a dictionary of those counts
in the order A, C, G, T or N if base is unknown

Examples  
≡≡≡≡≡≡≡≡≡≡

    julia> composition("AATCGGG")
    (2, 1, 3, 1)

    julia> composition('C')
    (0, 1, 0, 0)

    julia> A,C,G,T = composition("accgggtttt")
    (1, 2, 3, 4)

    julia> A
    1

    julia> T
    4

    julia> composition("BACCGGGTTTT")
    ERROR: Invalid base B encountered    
"""

sequence = "AATAACGGGNNTTN"
function composition(sequence)
    sequence = normalizeDNA(sequence) # make uppercase string, check invalid bases
    a = c = g = t = n = 0 # sets all 5 variables to `0`

    comp_dict = Dict()
    for base in sequence
        (base == 'A') && (a += 1)
        (base == 'T') && (t += 1)
        (base == 'C') && (c += 1)
        (base == 'G') && (g += 1)
        (base == 'N') && (n += 1)
        ## add 1 to each base as it occurs
    end
    comp_dict['A'] = a
    comp_dict['T'] = t
    comp_dict['C'] = c
    comp_dict['G'] = g
    comp_dict['N'] = n
    
    return comp_dict
end

function gc_content(sequence)
    gc = ["GC"...]
    gc_count = 0
    gc_content = float(0)
    sequence = uppercase(sequence)
    for letter in sequence
        !(letter in ["ATCGN"...]) && error("The sequence seems to include characters that are not DNA bases")
        (letter in gc) && (gc_count += 1)
    end
    (gc_count == 0) ? (gc_content = gc_count) : (gc_content = gc_count/length(sequence))
    return gc_content 
    ## Start with the same code as `question3()` from assignment 2.
    ## only a small modification is necessary to make this work.
end

function complement(seq)
    complements = Dict("A" => "T",
                       "T" => "A",
                       "G" => "C",
                       "C" => "G",
                       "N" => "N")

    complement_seq = ""
    for base in seq
        base = uppercase(string(base))
        !(base in keys(complements)) && error("Invalid base $base")
        complement_seq = complement_seq * complements[base]
    end
    return complement_seq
end

function reverse_complement(sequence)
    rev_seq_array = split(reverse(sequence),"")
    rev_comp_array = map(complement,rev_seq_array)
    rev_comp = join(rev_comp_array)
    return rev_comp
end

function parse_fasta(path)

    headers = []
    sequences = []
    subsequences = ""
    for line in eachline(path)
        ('>' in line) && push!(headers, string(strip(line, '>')))            
    end
    fasta_lines = readlines(path)
    for i in 1:(length(fasta_lines))
        if ('>' in fasta_lines[i] && subsequences != "")
            # println(subsequences)
            push!(sequences, normalizeDNA(subsequences))
            subsequences = ""
        elseif ('>' in fasta_lines[i] && subsequences == "")
            continue
        else
            subsequences = subsequences * fasta_lines[i]
        end
    end
    push!(sequences, normalizeDNA(subsequences))
    return (headers,sequences)
    ## Think through the components you need
    ## Does it make sense to define any containers at the beginning?
    ## How will you loop through the file?
    ## What do you need to get from each line?
end


# Don't forget to export your functions!


end # module Assignment07
