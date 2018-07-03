#!/usr/bin/env python
"""
Author: Stephen Coleman
Implementation of the SSAHA algorithm
"""
# import statements
from sys import argv
from itertools import product, chain, islice
from collections import Counter
from  operator import itemgetter
import re

BASES = "ACGT"

def remove_redundant_kmers(kmer_list, sequence, k):
    """Remove kmers not present in sequence from kmer_list

    kmer_list: list of strings (i.e. kemrs)
    sequence: str, sequence of DNA
    """
    final_list = set()
    for i in range(0, len(sequence)):
        if sequence[i : i + k] in kmer_list:
            final_list.add(sequence[i : i + k])
    return(list(final_list))

def list_kmers(sequence, k):
    """Returns all kmers contained in sequence

    sequence: str, sequence of DNA bases
    k: int, length of kmers
    """
    if not isinstance(k, int):
        raise(TypeError("k must be type int."))

    if not isinstance(sequence, str):
        raise(TypeError("sequence must be of type str"))

    kmer_list = set() # list()
    for i in range(0, len(sequence) - k + 1, k):
        kmer_list.add(sequence[i: i + k])

    if len(sequence) % k != 0:
        kmer_list.add(sequence[-k:])
    return(kmer_list)

def list_kmers_database(list_sequences, k):
    """Returns all kmers for a list of sequences.

    list_sequences: lst, database of sequences
    k: int, kmer length
    """
    kmer_list =set()
    for sequence in list_sequences:

        kmer_list = kmer_list.union(list_kmers(sequence, k))

    return(kmer_list)

def check_only_bases(sequence):
    """Tests if sequence contains any illegal symbols

    sequence: str, sequence of DNA bases
    """
    if not isinstance(sequence, str):
        raise(TypeError("sequence must be of type str"))


    set_bases = set(BASES)
    set_sequence = set(sequence)

    if not set_sequence.issubset(set_bases):
        return(False)
    return(True)

def remove_prohibited_symbols(sequence, allowed = "AGCT", replacement = "A"):
    """Checks for and removes all symbols not in allowed

    sequence: str, string being checked for prhibited symbols
    allowed: str, allowed characters within sequence (default "AGCT")
    replacement: str, character to replace forbidden characters with (default
    is 'A' in keeping with SSAHA paper)
    """
    if not isinstance(sequence, str):
        raise(TypeError("sequence must be of type str"))

    if not isinstance(allowed, str):
        raise(TypeError("allowed must be of type str"))

    if not (isinstance(replacement, str) and len(replacement) == 1):
        raise(TypeError("replacement must be a string of length 1"))

    regex = r"([^{}]+)".format(allowed)

    # Make sure sequence is only capital letters
    sequence = sequence.upper()

    while re.search(regex, sequence) is not None:
        match = re.search(regex, sequence)
        sequence = sequence[:match.start()] \
                    + replacement * (match.end() - match.start()) \
                    + sequence[match.end():]
    return(sequence)

def list_pointers(kmer, sequence):
    """Returns list of indices at which kmer begins in sequence

    kmer: str, kmer of bases
    sequence: str, DNA sequence
    """
    if not isinstance(kmer, str):
        raise(TypeError("kmer must be of type str"))

    if False:
        sequence = remove_prohibited_symbols(sequence, allowed = BASES)
    else:
        sequence = sequence.upper()

    i = 0
    pointers = []
    while sequence.find(kmer, i, len(sequence)) != -1:
        pointers += [sequence.find(kmer, i)]

        i = pointers[-1] + 1

    return(pointers)


def list_hits(query_sequence, list_sequences, k):
    """Returns a list of lists of hits for each kmer of a query_sequence 
    across a searched list of sequences

    query_sequence: str, sequence of DNA searching for hits
    list_sequences: lst, list of sequences used as database to search
    k: int, length of kmers
    """
    kmer_list = list_kmers(query_sequence, k)

    H = {}
    for kmer in kmer_list:
        offset = query_sequence.find(kmer)
        H_hit = []
        for i, sequence in enumerate(list_sequences):
            pointers = list_pointers(kmer, sequence)

            # Use + 1 to match the original paper
            H_hit += list(map(lambda x, y, z: (x + 1, y - z + 1, y + 1),
                              [i] * len(pointers),
                              pointers,
                              [offset] * len(pointers)
                              )
                          )
        H[kmer] = H_hit

    return(H)


def order_hits(hits):
    """Returns htis from query in order, sorted by index, then shift

    hits: list of hits from list_hits function, e.g. for a list of n hits
    hits = [(index_0, shift_0, offset_0), ..., (index_n, shift_n, offset_n)]
    """
    hits = [item for sublist in list(hits.values()) for item in sublist]
    ordered_hits = sorted(hits, key = itemgetter(0, 1))

    return(ordered_hits)

def find_longest_hit(ordered_hits):
    """Returns a tuple of the start position of the longest match and its 
    length

    ordered_hits: lst, list of hits (output from order_hits)
    """
    matches = Counter([(index, shift) for (index, shift, offset) in \
                      ordered_hits]
                      )


    # Want to return all best hits
    matches_dict = dict(matches)
    maxValue = max(matches_dict.values())

    # This is meant to allow multiple optimal alignments (single best is 
    # given by commented line)
    longest_hits = [key for key in matches_dict.keys() 
                    if matches_dict[key]==maxValue
                    ]
    return(longest_hits)

def find_alignment(query_sequence, list_sequences, k, include_reverse = False):
    """Returns two strings of the ungapped alignment from the query sequence
    with the most identical match from the sequences in list_sequences. End 
    gaps are used to match the length of the alignments

    query_sequence: str, DNA sequence searching for alignments
    list_sequences: list of strings, the database of sequences used in the 
    search
    k: int, the length of the kmers
    """
    if include_reverse:
        list_sequences += [sequence[::-1] for sequence in list_sequences]
    hits = list_hits(query_sequence, list_sequences, k)


    ordered_hits = order_hits(hits)
    longest_matches = find_longest_hit(ordered_hits)
    alignments = []

    
    hits_of_interest = []
    for longest_match in longest_matches:

        hits_from_match = [item for item in ordered_hits if item[0:2] 
                            == longest_match
                          ]

        hits_from_match = sorted(hits_from_match, key=lambda elem: elem[2])

        # Find which of the comparison sequences the best alignment is on
        sequence_num = longest_match[0] - 1 # -1 to return to python indexing
        comparison_seq = list_sequences[sequence_num]

        # Find the start of the alignment on the two sequences
        start_pos_seq = longest_match[1] - 1
        start_pos_query = hits_from_match[0][2] - hits_from_match[0][1]

        # Length of the alignment    
        match_len = hits_from_match[-1][2] + k

        # Number of end gaps to align sequences
        num_opening_gaps = start_pos_seq - start_pos_query
        num_closing_gaps = num_opening_gaps + len(query_sequence) \
                            - len(comparison_seq)

        connection = connection_maker(hits_from_match, sequence_num, k, comparison_seq)

        align1, align2, connection = construct_alignment(query_sequence, 
                                                         comparison_seq,
                                                         num_opening_gaps, 
                                                         num_closing_gaps,
                                                         connection
                                                         )

        alignments += [(align1, align2, connection)]
        hits_of_interest += [hits_from_match]

        # Ho ho 
        break

    return alignments, hits_of_interest, ordered_hits

def connection_maker(hits, sequence_num, k, sequence):
    """Creates the exact match comparison string

    hits: hits associated with best match
    """
    
    hits_on_sequence = [hit[2] -1 for hit in hits]
    
    connection = [" "] * hits_on_sequence[0] \
                 + ["|"] * (hits_on_sequence[-1] - hits_on_sequence[0] + k) \
                 + [" "] * (len(sequence) - hits_on_sequence[0])

    connection = "".join(connection)
    return(connection)

def construct_alignment(seq1, seq2, num_opening_gaps, num_closing_gaps, connection):
    """Returns two strings of the alignment of the two inputted sequences, 
    using the num_opening_gaps and num_closing_gaps to correct any 
    non-overlapping sections

    seq1: str, 
    seq2: str
    num_opening_gaps: int, number of opening gaps to place in front of seq1 
    (negative number implies in front of seq2)
    num_closing_gaps: int, as above, but appended to alignments and negative
    placed on seq1
    """

    # Construct alignment as two seperate strings
    align1 = seq1
    align2 = seq2

    if num_opening_gaps < 0:
        align2 = "-" * (-1 * num_opening_gaps) + align2
        connection = "-" * (-1 * num_opening_gaps) + connection
    else:
        align1 = "-" * num_opening_gaps + align1

    if num_closing_gaps < 0:
        align1 += "-" * (-1 * num_closing_gaps)
    else:
        align2 += "-" * num_closing_gaps
        connection += "-" * num_closing_gaps

    return align1, align2, connection

def fasta_to_dict(fasta_file_name):
    """Writes FASTA file to a dictionary with each key the name
    of the sequence in the FASTA file and the corresponding value
    the associated sequence converted to a single string.

    Key arguments:
    fasta_file_name: str, filename of FASTA file that contains the 
    sequences for analysis.
    """
    if type(fasta_file_name) != str:
        raise TypeError('Inputs must be strings.')

    try:
        fasta_dict = {}
        first = True
        with open(fasta_file_name) as fasta:
                for line in fasta:
                    if re.match(r'^>', line):
                        key = line.strip()[1:]
                        fasta_dict[key] = []
                        loc_line = ''

                    else:
                        fasta_dict[key] += [line.strip()]
    except IOError:
        raise IOError('File inaccesible')
    except FileNotFoundError:
        raise FileNotFoundError('File does not exist.')
    except Exception as err:
        raise Exception("{}. \nExiting function.".format(err))
    finally:
        return(fasta_dict)

def represent_kmer_as_int(kmer):
    """Returns the E(w) score mentioned in the paper

    kmer: str, kmer of bases
    """
    E_w = 0
    scores = dict(zip(BASES, list(range(4))))
    k = len(kmer)
    for i, base in enumerate(kmer):
        E_w += 4 ** (k - i - 1) * scores[base]
    return E_w

def find_alignments_multiple_queries(list_queries, list_sequences, k,
                                     n = 10,
                                     # printing = True,
                                     include_reverse = False
                                     ):
    """Prints alignments for list of query sequences within database of 
    sequences

    list_queries: lst, list of DNA sequences
    list_sequences: lst, list of DNA sequences comprising the database
    k: int, length of kmers to use in search
    """
    alignments, hits = [], []
    original_database_len = len(list_sequences)
    for i, query in enumerate(list_queries):
        print("Query sequence {}".format(i + 1))
        if include_reverse:
            list_sequences = list_sequences[:original_database_len] # as include reverse doubles this

        print_alignment(query, list(list_sequences.values()), k,
                        n = n,
                        printing_full = False,
                        include_reverse = include_reverse,
                        sequence_names = list(list_sequences.keys())
                        )

    return None, None

def relevant_alignment(connection, original_match_sequence, align2, k, n = 1):
    """Returns the section of alignment of interest

    """
    seq_pos = align2.find(original_match_sequence)
    start_pos = max(0, seq_pos + connection[0][2] - n)
    end_pos = min(seq_pos + connection[-1][2] + k + n, len(align2))

    return start_pos, end_pos


def print_alignment(query, list_sequences, k,
                    n = 10,
                    printing_full = True,
                    include_reverse = False,
                    sequence_names = None):

    alignments, hits_of_interest, ordered_hits = find_alignment(query,
                                                                list_sequences, 
                                                                k = k,
                                                                include_reverse = include_reverse
                                                                )

    for i, alignment in enumerate(alignments):
        connection = alignment[2]

        sequence_num = hits_of_interest[i][0][0]
        sequence_start = hits_of_interest[i][0][2]

        query_start = hits_of_interest[i][0][2] - hits_of_interest[i][0][1]
        
        sequence_end = hits_of_interest[i][-1][2] + k
        query_end = hits_of_interest[i][-1][2] - hits_of_interest[i][-1][1] + k

        start_pos, end_pos = relevant_alignment(hits_of_interest[i], 
            seqs[hits_of_interest[i][0][0] - 1] , alignment[1], k, n = n)

        print("\n")

        print("Number of hits across all sequences: {}".\
               format(len(ordered_hits))
               )

        hits_on_sequence = [hit for hit in ordered_hits if hit[0] == sequence_num + 1]
        print("Number of hits on matched sequence: {}".\
               format(len(hits_on_sequence))
               )

        if sequence_names is not None:
            search_name = sequence_names[sequence_num]
            print("Search sequence: {}".format(search_name))
        else:
            print("Search sequence number: {}".format(sequence_num))

        print("Match begins at base {:,} (query) and {:,} (search sequence)".\
               format(query_start, sequence_start))
        print("Match ends at base {:,} (query) and {:,} (search sequence)".\
               format(query_end, sequence_end))
        
        
        print("Query length: {:,}".format(len(query)))

        print("Alignment begins at hit: {}\nAlignment ends with hit: {}".\
               format(hits_of_interest[i][0], hits_of_interest[i][-1])
               )

        print("Alignment length: {:,}".format(hits_of_interest[i][-1][-1] 
                                              - hits_of_interest[i][0][-1] 
                                              + k
                                              )
              )

        if printing_full:
            for i in range(start_pos, end_pos, 80):
                print("{}\n{}\n{}\n\n".format(alignment[0][i : i + 80],
                                              connection[i : i + 80],
                                          alignment[1][i : i + 80]))

        print("\n")
        return None

def parsed_fasta_to_search_space(parsed_sequences):
    """Combines the parsed n sequences into n strings and counts total 
    number of bases and actual sequence

    parsed_sequences: dict of lists of strings, output of fasta_to_dict
    """
    number_of_bases_in_db = 0
    final_combined_seq = ""

    database = {}
    for sequence in parsed_sequences.keys():
        curr_seq_length = len("".join(parsed_sequences[sequence]))
        database[sequence] = "".join(parsed_sequences[sequence])
        print("Sequence name: {}".format(sequence))
        print("{} length: {:,}".format(sequence, curr_seq_length))
        number_of_bases_in_db += curr_seq_length
        final_combined_seq += "".join(parsed_sequences[sequence])

    return(database, number_of_bases_in_db, final_combined_seq)


def print_kmer_hash_table(sequences, k):
    """Prints hash table for given list of sequences and k

    sequences: list of strings
    k: int, length of kmers
    """
    hash_table_kmers = list_kmers_database(seqs, k)
    print("Kmer\tE(w)\t\t\t\t\tPositions")
    hash_table = {}
    scores = []
    for kmer in hash_table_kmers:
        hash_table[kmer] = []
        scores += [represent_kmer_as_int(kmer)]

    kmer_score_ordered = [kmer[0] for kmer in 
                          sorted(list(zip(hash_table_kmers, scores)), 
                                 key=lambda score: score[1])
                         ]

    for kmer in kmer_score_ordered:
        score = represent_kmer_as_int(kmer)
        for i, sequence in enumerate(seqs):
            pointers = list_pointers(kmer, sequence)
            pointers = [pointer + 1 for pointer in pointers]
            hits = list(zip([i + 1] * len(pointers), pointers))

            hash_table[kmer] += hits


        string_hits = ["({}, {})".format(hit[0], hit[1]) for hit in hash_table[kmer]]
        print("{}:\t{}\t{}".format(kmer, score, "\t".join(string_hits)))


if __name__ == "__main__":

    # the code below should produce the results necessary to answer 
    # the questions. In other words, if we run your code, we should 
    # see the data that you used to answer the questions
    
    # Set this to true to include reverse of database sequences in search space
    check_reverse = False

    query = 'TGCAACAT'
    
    s1 = 'GTGACGTCACTCTGAGGATCCCCTGGGTGTGG'
    s2 = 'GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT'
    s3 = 'GGATCCCCTGTCCTCTCTGTCACATA'
    seqs = [s1,s2,s3]

    k = 2

    print("\n=== Question 1 =================================")
    print("Hash table for k = {}".format(k))
    print_kmer_hash_table(seqs, k)

    print("\n=== Question 2 =================================")
    print_alignment(query, seqs, k, n = 30)


    print("\n=== Question 4 =================================")
    arabidopsis = fasta_to_dict(argv[1])

    database, arabidopsis_length, arabidopsis_seq =\
        parsed_fasta_to_search_space(arabidopsis)
    
    print("\nArabidopsis genome length: {:,}".format(arabidopsis_length))


    print("\n=== Quesiton 5 =================================")
    k = 15
    arabidopsis_hash_table = list_kmers_database(database.values(), k)
    print("Number of unique, non-overlapping kmers (k = {}) in arabadopsis: {:,}".\
          format(k, len(arabidopsis_hash_table))
          )


    print("\n=== Quesiton 6 =================================")

    queries = fasta_to_dict(argv[2])
    queries = ["".join(seq) for seq in queries.values()]
    
    if check_reverse:
        database += [sequence[::-1] for sequence in database]

    q1, hits_of_interest = find_alignments_multiple_queries(queries, database,
                                                            k = k,
                                                            include_reverse = False
                                                            )
