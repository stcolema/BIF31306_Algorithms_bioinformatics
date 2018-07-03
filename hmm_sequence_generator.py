#!/usr/bin/env python

"""
Author: Stephen D. Coleman (the D is for David)
Studnet nr: 940309160050
Dedication: This script goes out to all of the boys
Description: this is a script to ... fill the meaningless void :D
"""
# Import statements
from sys import argv
from random import random, randint
import re
from copy import deepcopy
from numpy.random import choice

# Function definitions

# Background amino acid probabilities
pa = { 'A':0.074, 'C':0.025, 'D':0.054, 'E':0.054, 'F':0.047, 'G':0.074,\
    'H':0.026, 'I':0.068, 'L':0.099, 'K':0.058, 'M':0.025, 'N':0.045,\
    'P':0.039, 'Q':0.034, 'R':0.052, 'S':0.057, 'T':0.051, 'V':0.073,\
    'W':0.013, 'Y':0.034 } 

# These are not summing to 1 due to rounding errors; for choice need
# the probabilities to sum closer to one, hence:
pa_weights = {key: prob / sum(pa.values()) for key, prob in pa.items()}

# List of possible Amino Acids (AAs)
AAs = list(pa.keys())

class HMM():
    """HMM object to store an HMM model

    This object is designed to keep track of all HMM states, emissions, and 
    transitions. It may be used in your implementation, but may also be 
    ignored, and replaced by a data structure of choice
    """
    # Emission probabilities for the match and insert states
    e_m   = []; e_i   = pa; 
    
    # Transition probabilities from/to matches, inserts and deletions
    t_mm  = []; t_mi  = []; t_md = [];
    t_im  = []; t_ii  = []
    t_dm  = []; t_dd  = []; 
    
    def __init__(self, nmatches):
        """Initialize HMM object with number of match states
        
        nmatches: int, number of match states
        """
    
        self.nmatches = nmatches
        
        self.e_m   = [dict(pa) for i in range(0,nmatches)]
        for i in range(0,nmatches):
            for j in pa.keys():
                self.e_m[i][j] = 0.0
        self.e_i   = pa;

        self.t_mm  = [0.0 for i in range(0,nmatches+1)]
        self.t_mi  = [0.0 for i in range(0,nmatches+1)]
        self.t_im  = [0.0 for i in range(0,nmatches+1)]
        self.t_ii  = [0.0 for i in range(0,nmatches+1)]
        self.t_md  = [0.0 for i in range(0,nmatches+1)]
        self.t_dm  = [0.0 for i in range(0,nmatches+1)]
        self.t_dd  = [0.0 for i in range(0,nmatches+1)]

def sample(events):
    """Return a key from dict based on the probabilities 

    events: dict of {key: probability}, sum of probabilities should be 1.0. 
    """
    key_options = list(events.keys())
    cum = [0 for i in key_options]

    cum[0] = events[key_options[0]]
    for i in range(1, len(events)):
        cum[i] = cum[i-1] + events[key_options[i]]

    # Should not be necessary, but for safety
    cum[len(cum)-1] = 1.0

    ref_point = random()
    pick = None
    i = 0
    while not pick and (i < len(cum)):
        if ref_point < cum[i]:
            pick = key_options[i]
        i = i + 1
    return pick
    
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

def collapse_fasta_dict_entries_to_strings(fasta_dict):
    """Returns dict with values changed from lists of strings to one string

    fasta_dict: dict of lists of strings (output of fasta_to_dict)
    """
    for key in fasta_dict.keys():
        fasta_dict[key] = "".join(fasta_dict[key])
    return fasta_dict

def most_common(lst):
    """Returns the most common element in a list with its count

    lst: a list
    """
    return max(set(lst), key=lst.count)

def common_sequence_length(sequences):
    """Check if all sequences of common length

    sequences: list of strings representing different sequences
    """
    lengths = set([len(sequence) for sequence in sequences])
    if len(lengths) != 1:
        return False
    return True

def estimate_match_states(sequences, threshold = 0.5):
    """Returns the approximate number of match states in a list of sequences

    sequences: list of strings
    threshold: float; the fraction of the number of sequences that must match
    before the index in the sequences is considered a possible match state
    (default = 0.5)
    """
    if not common_sequence_length(sequences):
        raise ValueError("Not all sequences of common length.")

    # Initialise the number of match states
    num_match_states = 0

    # Cycle through index of sequences checking for a match state at each 
    # point
    for i in range(len(sequences[0])):

        # Check if a match state in the current position
        if check_match_state(sequences, i, threshold = threshold):
            num_match_states += 1

    return num_match_states

def residues_present(sequences, index):
    """Returns a list of the AAs present at current index in sequences

    sequences: list of strings
    index: int; index of strings in sequences to check as a match state
    """
    # Finds the residues in the current position across each sequence
    aas_present = [sequence[index] for sequence in sequences 
                   if sequence[index].upper() in AAs
                   ]

    return aas_present

def check_match_state(sequences, index, threshold = 0.5):
    """Checks if the current index in the list of sequences is a match state

    sequences: list of strings
    index: int; index of strings in sequences to check as a match state
    threshold: float; the fraction of the number of sequences that must match
    before the index in the sequences is considered a possible match state
    (default = 0.5)
    """
    # Finds the residues in the current position across each sequence
    aas_present = [sequence[index] for sequence in sequences 
                   if sequence[index].upper() in AAs
                   ]

    # The number of AAs present in the current index
    frequency = len(aas_present)
    
    # If this exceeds our threshold for defining a match state, add 1 to
    # the count

    if frequency >= threshold * len(sequences):
        return True
    return False

def current_states(sequences, 
                   threshold = 0.5,
                   match = "match",
                   insertion = "insertion",
                   deletion = "deletion"
                   ):
    """Returns the state for each sequence in sequences

    sequences: list of strings representing AA sequences
    index: int; index of strings in sequences to check as a match state
    threshold: float; the fraction of the number of sequences that must match
    before the index in the sequences is considered a possible match state
    (default = 0.5)
    match: str; tag identifying matches in path
    insertion: str; tag identifying insertions in path
    deletion: str; tag identifying deletions in path
    """

    states = []
    match_state = []

    for i, query in enumerate(sequences):
        states += [[]]

        match_state_num = 0

        for index in range(len(query)):

            # We do this seperate method for the 0th element to avoid 
            # repeating the check_match_state()
            if i == 0:

                # If a match state, we can be an insertion or a match
                if check_match_state(sequences, index, threshold):

                    # This number is used to update our current step, so if
                    # no match states are hit, each insertion will record
                    # the same number (i.e. we remain at insertion_i until
                    # we reach the (i + 1)th match state)
                    match_state_num += 1

                    match_state += [True]
                    
                    # If an amino acird is present, record a match
                    if query[index] in AAs:
                        states[i] += [match + "_{}".format(match_state_num)]
                        continue

                    # Else record a deletion
                    else:
                        states[i] += [deletion + "_{}".format(match_state_num)]
                        continue

                # If not a match state can be a insertion only (or else a gap 
                # present which is not recorded)
                else:
                    match_state += [False]

                    if query[index] in AAs:

                        states[i] += [insertion + "_{}".format(match_state_num)]

            # Logic as for i==0 but we have recorded if current step is a
            # match state or not
            else:
                if match_state[index]:
                    match_state_num += 1
                    if query[index] in AAs:
                        states[i] += [match + "_{}".format(match_state_num)]
                        continue

                    else:
                        states[i] += [deletion + "_{}".format(match_state_num)]
                        continue

                if query[index] in AAs:
                    states[i] += [insertion + "_{}".format(match_state_num)]
    return states

def count_transitions(sequence_states, insertion = "insertion"):
    """Returns the number of transitions between the states of a sequence

    sequence_states: list of ordered ints representing states
    insertion: str; tag identifying insertions in path
    """
    transitions = 1 # Account for start, should we have one for end?
    prev = sequence_states[0]

    for state in sequence_states[1:]:
        if not (state == prev and insertion in state):
            transitions += 1
        prev = state
    return transitions

def number_transitions(sequences, threshold = 0.5):
    """Returns a list of counts for the number of transitions in each sequence

    sequences: list of strings representing AA sequences
    threshold: float; fraction of sequences which must have a AA at a given
    index for said index to be considered a match state
    """
    total_transitions = []
    states = current_states(sequences, threshold = threshold)
    # for state in states:
    #     print(convert_state_int_to_string(state))
    for sequence_states in states:
        total_transitions += [count_transitions(sequence_states)]
    return total_transitions

def initialise_states(sequences,
                      match = "match",
                      insertion = "insertion",
                      deletion = "deletion"
                      ):
    """Initialises the possible states available for a given set of sequences

    sequences: list of strings representing AA sequences
    match: str; tag identifying matches in path
    insertion: str; tag identifying insertions in path
    deletion: str; tag identifying deletions in path
    """
    # Find the number of match states
    num_match_states = estimate_match_states(sequences)

    # Declare the possible match, insertion and deletion states
    matches = ["{}_{}".format(match, i) 
               for i in range(1, num_match_states + 1)
               ]

    deletes = ["{}_{}".format(deletion, i) 
               for i in range(1, num_match_states + 1)
               ]

    # Note the different range for insertions
    inserts = ["{}_{}".format(insertion, i) 
               for i in range(0, num_match_states + 1)
               ]

    # Start node
    start = ["start"]

    possible_states = start + matches + inserts + deletes #+ end

    # Initial entries are empty lists
    entries = [[] for i in possible_states]

    # Create a dict (treating as a graph)
    initial_map = dict(zip(possible_states, entries))

    return initial_map

def transition_path(sequences,
                    match = "match",
                    insertion = "insertion",
                    deletion = "deletion"
                    ):
    """Returns the path through transitions for a set of sequences

    sequences: list of strings representing AA sequences
    match: str; tag identifying matches in path
    insertion: str; tag identifying insertions in path
    deletion: str; tag identifying deletions in path
    """
    inital_states = initialise_states(sequences, match, insertion, deletion)
    states = current_states(sequences)
    for state in states:
        inital_states["start"] += [state[0]]
        for i, transition in enumerate(state):
            if i < len(state) - 1:
                inital_states[transition] += [state[i + 1]]
            else:
                inital_states[transition] += ["end"]
    return inital_states

def weighted_path(sequences,
                  match = "match",
                  insertion = "insertion",
                  deletion = "deletion"
                  ):
    """Returns a dict of dicts key corresponds to state, entries to transitions

    sequences: list of strings representing AA sequences
    match: str; tag identifying matches in path
    insertion: str; tag identifying insertions in path
    deletion: str; tag identifying deletions in path
    """
    transitions = transition_path(sequences,
                                  match = match,
                                  insertion = insertion,
                                  deletion = deletion
                                  )

    normalised_path = deepcopy(transitions)
    for state in transitions.keys():
        weights = []
        options = set(transitions[state])
        for option in options:
            weights += [transitions[state].count(option)
                        / len(transitions[state])
                        ]
        curr_entry = dict(zip(options, weights))
        normalised_path[state] = curr_entry
    return normalised_path

def realised_emissions(sequences):
    """Returns a list of the observed emissions for each match state

    sequences: list of strings representing AA sequences
    """
    emissions = []
    # for sequence in sequences:
    for i, res in enumerate(sequences[0]):
        if check_match_state(sequences, i):
            emissions += [residues_present(sequences, i)]

    return emissions

def weigh_emissions(sequences, match = "match", pseudocounts = True):
    """Returns the weigths for each emission for each match state

    sequences: list of strings representing AA sequences
    match: str; tag identifying matches in path
    pseudocounts: bool; instructs usage of pseudocounts in match emission
    """
    emissions = realised_emissions(sequences)
    emission_weights = {}

    for i, elem in enumerate(emissions):
        weights = []
        if pseudocounts:
            for aa in AAs:
                weights += [(elem.count(aa) + 1)/ (len(elem) + 20)]

            curr_entry = dict(zip(AAs, weights))

        else:
            options = set(elem)
            for option in options:
                weights += [elem.count(option) / len(elem)]

            curr_entry = dict(zip(options, weights))
        emission_weights["{}_{}".format(match, i + 1)] = curr_entry
    return emission_weights

def generate_sequence(sequences, 
                      match = "match", 
                      insertion = "insertion",
                      deletion = "deletion",
                      pseudocounts = True,
                      silent_state = True
                      ):
    """Generates a sequence and path from sequences based on HMM

    match: str; tag identifying matches in path
    insertion: str; tag identifying insertions in path
    deletion: str; tag identifying deletions in path
    pseudocounts: bool; instructs usage of pseudocounts in match emission
    """

    # Generate HMM based on sequences; here we identify possible paths and 
    # emissions
    sequence_path = weighted_path(sequences, match, insertion, deletion)
    match_emissions = weigh_emissions(sequences,
                                      match = match,
                                      pseudocounts = pseudocounts
                                      )
    insert_emissions = pa_weights

    # Initialise the output path and sequence
    new_sequence_path = []
    new_sequence = []

    # Step through path, recording emissions and steps in each iteration
    step = 'start'
    while step != 'end':
        options = list(sequence_path[step].keys())
        probabilities = list(sequence_path[step].values())
        next_step = choice(options, p = probabilities)
        if match in next_step:

            # Choose from options using weighted probabilities
            emission = choice(list(match_emissions[next_step].keys()),
                              p = list(match_emissions[next_step].values())
                              )
        elif insertion in next_step:
            emission = choice(list(insert_emissions.keys()),
                              p = list(insert_emissions.values())
                              )
        elif deletion in next_step and not silent_state:
            emission = "-"
        if next_step != "end":
            new_sequence += [emission]
        new_sequence_path += [step]
        step = next_step

    # Record end step
    new_sequence_path += [step]

    # Convert sequence from a list of characters to a string
    new_sequence = "".join(new_sequence)

    # Return the path that generated the sequence and the sequence
    return new_sequence_path, new_sequence

def print_PSSM_matrix(emission_weights):
    """Print a matrix of position specific probabilities for each AA

    emission_weights: dict of dicts of the possible emissions and their
    associated probabilities for each match state
    """

    # keep output as a string and then join before printing
    output = [' ', ' ', ' '] + AAs + ['\n']
    for i, match in enumerate(emission_weights.keys()):
        new_line = [str(i), "A"]
        for aa in AAs:
            if aa in emission_weights[match].keys():
                curr_prob = emission_weights[match][aa]
            else:
                curr_prob = 0.0

            new_line += ["{:.3f}".format((curr_prob))]
        new_line += ['\n']
        output += new_line

    out = "\t".join(output)
    print(out)

def requested_output(infile, num_generated_sequences = 10, save_seq = False):
    """Prints the outputs requested in assignment

    infile: str; filename containing aligned sequences to parametrise a HMM
    num_generated_sequences: int; number of sequences to generate from HMM
    (default is 10)
    save_seq: bool; instructs printing of a random sequence from generated
    sequences, stripped of gaps (i.e. ready for use in PFAM search)
    """
    # Read in the fasta file
    sequences = fasta_to_dict(infile)

    # Collapse entries to single strings rather than lists of strings
    sequences = collapse_fasta_dict_entries_to_strings(sequences)
    seqs = list(sequences.values())

    # Print a seperating line naming file (length is to match PSSM matrix)
    additional_breaks = "=" * max((169 - len(infile)),0)
    print("\n === File: {} {}".format(infile, additional_breaks))
    
    # Estimate the number of match states
    print("\nNumber of match states: {}".\
            format(estimate_match_states(seqs))
          )


    # If want to save a sequence to search PFAM database, randomly select
    # using (0, 1) for whether in pseudocount iteration or not, and then the
    # sequence number
    if save_seq:
        choice = (randint(0, 1), randint(0, num_generated_sequences - 1))

    # Iterate through use of pseudocounts
    for i, pseudocounts in enumerate([False, True]):

        print("\nUsing pseudocounts: {}".format(pseudocounts))
        
        emission_weights = weigh_emissions(seqs, pseudocounts = pseudocounts)
        transition_path = weighted_path(seqs,
                                        match = "match",
                                        insertion = "insertion",
                                        deletion = "deletion"
                                        )
    
        print_transition_path(transition_path)
        print_emission_weights(emission_weights)

        # Print PSSM matrix
        print("\nPSSM matrix:")
        print_PSSM_matrix(emission_weights)
        
        # Print the sequences generated using the HMM
        print("Generated sequences:")
        for j in range(num_generated_sequences):

            gen_seq = generate_sequence(seqs, pseudocounts = pseudocounts)
            if save_seq:
                if i == choice[0] and j == choice[1]:
                    saved_sequence = gen_seq[1]
            print(gen_seq[1])

    if save_seq:
        print("\nSaved sequence:")
        saved_sequence = saved_sequence.replace("-", "")
        print(saved_sequence)

def print_transition_path(transition_path):
    """Prints the possible paths and associated probabilities through a HMM

    transition_path: dict of dicts containing a graph with weighted edges
    key (level 1) string describing current state (e.g. "start", "match_12"),
    value (level 1) dict describing possible movements from current state
    key (level 2) target state as a str and value (level 2) associated 
    probability of moving to this state
    """
    for step in transition_path.keys():
        if transition_path[step]:
            print("\n{}".format(step.upper()))
            out_str = [' '] + list(transition_path[step].keys()) + ['\n'] \
                      + ["{:.3}".format(prob) 
                         for prob in transition_path[step].values()
                         ]
            print("\t".join(out_str))

def print_emission_weights(emission_weights):
    """Prints possible emission and associated probabilities for match states

    emission: dict of dicts containing a graph with weighted edges describing
    key (level 1) the match state, value (level 1) dict
    key (level 2) the emitted AA, value (level 2) probability as a float
    """
    for match in emission_weights.keys():
        print("\n{}".format(match.upper()))
        out_str = [' '] + list(emission_weights[match].keys()) + ['\n'] \
                  + ["{:.3f}".format(prob) 
                     for prob in emission_weights[match].values()
                     ]

        print("\t".join(out_str))

if __name__ == "__main__":

    # implement main code here
    infile = 'test.fasta'
    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file) 
    requested_output(infile)

    # Second file and output
    infile = 'test_large.fasta'
    requested_output(infile, save_seq = True)