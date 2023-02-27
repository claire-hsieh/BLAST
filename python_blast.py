# python version of simple BLAST alg

query = "didyoseeasteelrvlashorwasitaflagthere"
seq = "mkvlwaallvtflagcqakveqavetepepelrqqtewqsgqrwelalgrfwdylrwvqtlseqvqeellssqvtqelralmdetmkelkaykseleeqltpvaeetrarlskelqaaqarlgadvlashgrlvqyrgevqamlgqsteelrvrlashlrklrkrllrvlashqkrlavyqagaregaerglsairerlgplveqgrvraatvgslagqplqeraqawgerlrarmeemgsrtrdrldevkeqvaevrakleeqaqqirlvlashqarlkswfeplvedmqrqwaglvek"
k = 4

# get substrings of query
queries = []
for i in range(len(query)-k+1):
    queries.append(query[i:(i+k)])

# get substrings of the sequence
sequences = []
for j in range(len(seq)-k+1):
    sequences.append(seq[j:(j+k)])

# split query / seq into list of characters
q_char = list(query)
s_char = list(seq)

# create database of location of kmers
s_database = {}
x = []
for ind, s in enumerate(sequences):
    if s in s_database:
        x = s_database[s][0:]
        x.append(ind)
        s_database[s] = x
        x = []
    else:
        s_database[s] = [ind]
        
string_hash = {}
# iterate through kmers in query and check if in database
for q_index, l in enumerate(queries):
    l_length = 0
    ls_l_lengths = []
    r_length = 0
    ls_r_lengths = []

    if l in s_database:
        if len(s_database[l]) > 1:
            # iterate thru all loc of kmer in seq
            # use loc with longest match
            # check left
            for s_index in s_database[l]:
                while (q_index >= 0 and s_index >= 0):
                    if q_char[q_index] == s_char[s_index]:
                        l_length += 1
                        q_index -= 1
                        s_index -= 1
                    else:
                        ls_l_lengths.append((l_length, s_index+l_length))
                        l_length = 0
                        q_index = queries.index(l)
                        break
            # check right
            s_index = sequences.index(l) + k
            q_index = queries.index(l) + k
            for s_index in s_database[l]:
                while (q_index <= len(queries) and s_index <= len(sequences)):
                    if q_char[q_index] == s_char[s_index]:
                        r_length += 1
                        q_index += 1
                        s_index += 1
                    else:
                        ls_r_lengths.append((r_length, s_index-r_length+1))
                        r_length = 0
                        q_index = queries.index(l)+ k
                        break
            # get longest match
            s_index = 0
            max_match = 0
            for i, j in zip(ls_l_lengths, ls_r_lengths):
                if (i[0] + j[0]) > max_match:
                    l_length = i[0]
                    r_length = j[0]
                    s_index = i[1]
                    max_match = l_length + r_length
        else:
            s_index = s_database[l][0]
            # check left
            while (q_index >= 0 and s_index >= 0):
                if q_char[q_index] == s_char[s_index]:
                    l_length += 1
                    q_index -= 1
                    s_index -= 1
                else:                
                    break

            # check right
            s_index = s_database[l][0] + k
            q_index = queries.index(l) + k
            while (q_index <= len(queries) and s_index <= len(sequences)):
                if q_char[q_index] == s_char[s_index]:
                    r_length += 1
                    q_index += 1
                    s_index += 1
                else:
                    break
        length = l_length + k + r_length
        if length >= k:
            output = "".join(s_char[sequences.index(l)-l_length+1: sequences.index(l)+k+r_length])
            if output not in string_hash:
                string_hash[output] = 1
                print(l, output, length, sequences.index(l))

            