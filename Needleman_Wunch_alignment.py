#!/usr/bin/env python3

"""
Author: Stephen Coleman

Description: this is a script to implement the Needleman Wunsch algorithm for 
aligining two sequences of peptides using a BLOSUM62 scoring matrix.
"""
#import statements here
from numpy import argmax
# import sys # no longer used

# functions between here and __main__
blosum = """
# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""

def blosum62():
    """Return order and similarity scores from BLOSUM62 matrix

    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {}
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 24:
            for idx, sym in enumerate(parts):
                order[sym] = idx
        else:
            # list around the map construction for python3 compatibility
            blosum_matrix.append(list(map(int,parts[1:])))
    return order, blosum_matrix

BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()

def score(res1, res2):
    """Return similarity score from BLOSUM62 matrix for two residues
    
    res1: string, amino acid
    res2: string, amino acid
    """
    lookup1 = BLOSUM62_ORDER[res1]
    lookup2 = BLOSUM62_ORDER[res2]
    return BLOSUM62_MATRIX[lookup1][lookup2]

def gap_penalty(number_gaps, initial, extension = None):
    """Return gap penalty; default is a linear gap penalty.

    number_gaps: int, number of gaps
    initial: penalty for initialising a gap
    extension: penalty for extending a gap (default is initial)
    """
    if not isinstance(number_gaps, int):
        raise(TypeError("number_gaps must be type int."))

    if not isinstance(initial, int):
        raise(TypeError("initial must be type int."))

    if extension == None:
        extension = initial

    if not isinstance(extension, int):
        raise(TypeError("extension must be type int."))

    return initial + (number_gaps - 1) * extension


def initialise_matrix(num_rows, num_cols, 
    initial_gap_penalty,
    extension_gap_penalty = None,
    end_penalty = None,
    semiglobal = False
    ):
    """ Initialises matrix for Needleman-Wunsch.

    num_rows: int, number of rows in matrix
    num_cols: int, number of columns in matrix
    initial_gap_penalty: int, initial gap penalty
    extension_gap_penalty: int, extension gap penalty
    end_penalty: bool, indicates a penalty for gaps at either end of alignment
    """
    if not isinstance(num_rows, int):
        raise(TypeError("num_rows must be type int."))

    if not isinstance(num_cols, int):
        raise(TypeError("num_cols must be type int."))



    if end_penalty is None:
        end_penalty = initial_gap_penalty

    if not isinstance(end_penalty, int):
        raise(TypeError("end_penalty must be type int."))   

    elif end_penalty != initial_gap_penalty:
        extension_gap_penalty = end_penalty
    matrix = [[0 for j in range(num_cols + 1)] for i in range(num_rows + 1)]
    if end_penalty != 0 and not semiglobal:
        for i in range(1, num_rows + 1):
            matrix[i][0] = - gap_penalty(i, 
                                         end_penalty,
                                         extension_gap_penalty
                                        )
        for j in range(1, num_cols + 1):
            matrix[0][j] = - gap_penalty(j, 
                                         end_penalty,
                                         extension_gap_penalty
                                        )

    return matrix

def initialise_path(num_rows, num_cols):
    """Initialise the array recording the path from beginning to end.

    num_rows: int, the number of rows in the alignment matrix
    num_cols: int, the number of columns in the alignment matrix
    """
    path = [[[0, 0] for j in range(num_cols + 1)] for i in range(num_rows + 1)]
    for i in range(1, num_rows + 1):
        path[i][0] = [-1, 0]
    for j in range(1, num_cols + 1):
        path[0][j] = [0, -1]
    return path


def Needleman_Wunsch(res1, res2, 
    initial_gap_penalty,
    extension_gap_penalty = None,
    end_penalty = None,
    semiglobal = False):
    """Needleman-Wunsch algorithm for comparing sequences using BLOSUM62

    res1: str, first sequence for alignment
    res2: str, second sequence for alignment
    initial_gap_penalty: int, penalty for beginning a gap
    extension_gap_penalty: int, penalty for continuing a gap 
    (defaults to the same as initial)
    end_penalty: penalty for having a gap at either end of the sequences
    (defaults to the same as initial)
    """
    extension_gap_penalty = None # as this flexibility isn't built in

    num_rows, num_cols = len(res1), len(res2)
    path = initialise_path(num_rows, num_cols)

    multiple_paths = False

    matrix = initialise_matrix(num_rows, num_cols, 
                               initial_gap_penalty = initial_gap_penalty,
                               extension_gap_penalty = extension_gap_penalty,
                               end_penalty = end_penalty,
                               semiglobal = semiglobal
                               )
    if end_penalty is None:
        end_penalty = initial_gap_penalty


    for i in range(1, num_rows + 1):
        for j in range(1, num_cols + 1):

            options = [matrix[i - 1][j - 1] 
                       + score(res1[i - 1], res2[j - 1]),
                       matrix[i - 1][j] - initial_gap_penalty,
                       matrix[i][j - 1] - initial_gap_penalty
                      ]

            # Account for end gap penalty in final row and column
            if i == num_rows or j == num_cols:
                options = [matrix[i - 1][j - 1] 
                           + score(res1[i - 1], res2[j - 1]),
                           matrix[i - 1][j] - end_penalty,
                           matrix[i][j - 1] - end_penalty
                          ]

            matrix[i][j] = options[argmax(options)]

            path[i][j] = [-1, -1] * int(argmax(options) == 0) \
                         + [-1, 0] * int(argmax(options) == 1) \
                         + [0, -1] * int(argmax(options) == 2)

    return matrix, path

def traceback(path, matrix = None, semiglobal = False):
    """Returns the path of alignment, favouring the diagonal.

    path: matrix of movements from NW function
    matrix: aignment matrix
    semiglobal: bool, True ignores end gap penalties
    """
    if semiglobal and matrix is None:
        raise(ValueError("matrix must be given if semiglobal."))

    num_rows = len(path)
    num_cols = len(path[0])

    if semiglobal:        
        max_j = max(matrix[num_rows - 1])

        final_col = [matrix[i][num_cols - 1] for i in range(num_rows)]
        max_i = max(final_col)

        if max_i < max_j:
            i, j = num_rows - 1, argmax(final_col)
        else:
            i, j = argmax(matrix[num_rows - 1]), num_cols - 1

    else:
        i, j = num_rows - 1, num_cols - 1

    
    final_path = [(i, j),]

    while ((i > 0 or j > 0) and not semiglobal) \
        or ((i > 0 and j > 0) and semiglobal):

        step = (path[i][j][0] + final_path[0][0],
                path[i][j][1] + final_path[0][1]
                )

        final_path = [step] + final_path
        i, j = step[0], step[1]

    return final_path

def print_NW_output(matrix, ref1, ref2):
    """Prints the alignment matrix with the sequences as borders

    matrix: alignment matrix
    ref1: str, sequence being compared in matrix
    ref2: str, second sequence compared in matrix
    """
    num_rows = len(matrix)
    num_cols = len(matrix[0])

    printed_matrix = [["" for j in range(num_cols + 1)] 
                      for i in range(num_rows + 1)
                     ]

    for i in range(2, num_rows + 1):
        printed_matrix[i][0] = ref1[i - 2]

    for j in range(2, num_cols + 1):
        printed_matrix[0][j] = ref2[j - 2]

    for i in range(0, num_rows):
        for j in range(0, num_cols):
            printed_matrix[i + 1][j + 1] = matrix[i][j]

    for row in printed_matrix:
        for elem in row:
            print(elem, end='\t')
        print()

def trace_alignment(seq1, seq2, trace):
    """ Returns two strings, the sequences including the gaps due to alignment.

    seq1: str, sequence being compared
    seq2: str, second sequence being compared
    trace: output of traceback function for seq1 and seq2
    """
    align1 = ""
    align2 = ""
    seq1_ind = trace[0][0]
    seq2_ind = trace[0][1]
    for i in range(1, len(trace)):
        if trace[i][0] > trace[i - 1][0] and trace[i][1] > trace[i - 1][1]:
            align1 += seq1[seq1_ind]
            align2 += seq2[seq2_ind]
            seq1_ind += 1
            seq2_ind += 1
        elif trace[i][0] > trace[i - 1][0] and trace[i][1] == trace[i - 1][1]:
            align1 += seq1[seq1_ind]
            align2 += "-"
            seq1_ind += 1
        else:
            align1 += "-"
            align2 += seq2[seq2_ind]
            seq2_ind += 1
    return align1, align2

def identity(align1, align2):
    """Calculates the identity between two alignments

    align1: str, sequence including gaps
    align2: str, second sequence being compared, including gaps
    """
    num_common = len(list(filter(lambda x: (x[0] == x[1] and x[0] != "-"), 
                                 zip(align1, align2)
                                )
                         )
                    )

    align_len = len(align1)
    return (num_common / align_len) * 100


def question_1_print(seq1, seq2, initial_gap_penalty, 
                     end_gap_penalty = None, 
                     semiglobal = False, 
                     max_print_len = 10
                    ):
    """Function to print output required by question 1.

    seq1: str, sequence to be compared
    seq2: str, sequence to be compared
    initial_gap_penalty: int, penalty for initialising a gap
    end_gap_penalty: int, penalty for a gap at the end of the alignment, 
    defaults to the same value as initial_gap_penalty
    semiglobal: bool, allows semiglobal alignment, default is False
    max_print_len: int, restricts number of characters from each sequence to 
    print
    """

    if end_gap_penalty is None:
        end_gap_penalty = initial_gap_penalty

    # Bool indicating if sequences exceed printable length
    too_long = False

    alignment_matrix, path = Needleman_Wunsch(seq1,seq2, 
                               initial_gap_penalty = initial_gap_penalty, 
                               end_penalty = end_gap_penalty,
                               semiglobal = semiglobal
                               )
    
    trace = traceback(path, alignment_matrix, semiglobal = semiglobal)

    align_seq1, align_seq2 = trace_alignment(seq1, seq2, trace)

    seq1_to_print = min(len(seq1), max_print_len)
    seq2_to_print = min(len(seq2), max_print_len)

    if max_print_len is None:
        pass
    elif max_print_len < len(seq1) or max_print_len < len(seq2):
        print("Abbreviating printed statements due to sequence length.")
        too_long = True

    print("\nScore matrix used: BLOSUM62")
    print("Initial gap penalty: {}\nEnd gap penalty: {}".\
        format(initial_gap_penalty, end_gap_penalty))

    print("Sequences being compared: \nSEQ1: {}\nSEQ2: {}".\
              format(seq1[:seq1_to_print], seq2[:seq2_to_print]))
    if not too_long:
        print_NW_output(alignment_matrix, seq1, seq2)
        print("Sequence alignment:\n{}\n{}\n".format(align_seq1, align_seq2))
        print("Path through matrix is:\n{}\n".format(trace))
    
    print("Score: {}".format(alignment_matrix[len(alignment_matrix) - 1]\
          [len(alignment_matrix[0]) - 1]))

    print("Identity: {0:.1f}%".format(identity(align_seq1, align_seq2)))
    
    print("#================================================================")

    return 0


def question_3(seq1, seq2, min_pen = 1, max_pen = 20, full_output = False):
    """Function to produce output relevant to Q3 of assignment. Finds all 
    unique alignments of seq1 and seq2 using a BLOSUM62 scoring matrix and 
    a gap penalty in the range(min_pen, max_pen + 1)

    seq1: str, first sequence being compared
    seq2: str, second sequence being compared
    min_pen: int, initial gap penalty used in iteration (default = 1)
    max_pen: int, final gap penalty used in iteration (default = 20)
    """
    pen = []
    original_traces = []
    original_alignments = []
    associated_matrices = []
    for initial_gap_pen in range(min_pen, max_pen + 1):
        alignment_matrix, path = Needleman_Wunsch(seq1,seq2, 
                               initial_gap_penalty = initial_gap_pen, 
                               end_penalty = initial_gap_pen
                               )
        current_traceback = traceback(path)
        if current_traceback not in original_traces:
            original_traces += [current_traceback]
            associated_matrices += [alignment_matrix]
            pen += [initial_gap_pen]
            align_seq1, align_seq2 = trace_alignment(seq1, 
                                                     seq2,
                                                     current_traceback
                                                     )

            original_alignments += [(align_seq1, align_seq2)]

    if full_output:
        for i in range(len(original_alignments)):

            question_1_print(seq1, seq2, pen[i], 
                             end_gap_penalty = None,
                             semiglobal = False
                             )


    print("\nUnique alignments for gap penalties [1, 20] and sequences:" \
          + "\n{}\n{}".format(seq1,seq2)
         )
    for i in range(len(original_alignments)):
        print("\n")
        print("Gap penalty: {}".format(pen[i]))
        print("Alignment: \n{}\n{}".format(original_alignments[i][0],
                                           original_alignments[i][1]
                                           )
              )
    print("\n")

    print("#================================================================")
    return 0

def question_4(seq1, seq2, initial_gap_penalty, end_gap_penalty):
    """Specific function qwustion 4 output, using question 1 function.
    Purely for ease of seeing where question 4 is done in code.
    """
    question_1_print(seq1, seq2, initial_gap_penalty, 
                     end_gap_penalty = end_gap_penalty
                     )
    return 0


if __name__ == "__main__":

    seq1 = "THISLINE"
    seq2 = "ISALIGNED"
    # seq3: GPA1_ARATH
    seq3 = "MGLLCSRSRHHTEDTDENTQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQIKLLFQTGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDSAKYMLSSESIAIGEKLSEIGGRLDYPRLTKDIAEGIETLWKDPAIQETCARGNELQVPDCTKYLMENLKRLSDINYIPTKEDVLYARVRTTGVVEIQFSPVGENKKSGEVYRLFDVGGQRNERRKWIHLFEGVTAVIFCAAISEYDQTLFEDEQKNRMMETKELFDWVLKQPCFEKTSFMLFLNKFDIFEKKVLDVPLNVCEWFRDYQPVSSGKQEIEHAYEFVKKKFEELYYQNTAPDRVDRVFKIYRTTALDQKLVKKTFKLVDETLRRRNLLEA"
    # seq4: GPA1_ORYSI
    seq4 = "MGSSCSRSHSLSEAETTKNAKSADIDRRILQETKAEQHIHKLLLLGAGESGKSTIFKQIKLLFQTGFDEAELRSYTSVIHANVYQTIKILYEGAKELSQVESDSSKYVISPDNQEIGEKLSDIDGRLDYPLLNKELVLDVKRLWQDPAIQETYLRGSILQLPDCAQYFMENLDRLAEAGYVPTKEDVLYARVRTNGVVQIQFSPVGENKRGGEVYRLYDVGGQRNERRKWIHLFEGVNAVIFCAAISEYDQMLFEDETKNRMMETKELFDWVLKQRCFEKTSFILFLNKFDIFEKKIQKVPLSVCEWFKDYQPIAPGKQEVEHAYEFVKKKFEELYFQSSKPDRVDRVFKIYRTTALDQKLVKKTFKLIDESMRRSREGT"
    
    print("## ====================== QUESTION 1 ============================")
    # OK so setting semiglobal to True gives odd results for identity score.
    # Let's not use semiglobal. But I am not sure how to approach this, would
    # comparing the full sequences be correct (as is in the case of 1c or any
    # 0 end penalty scenario) or should we compare the subsequences actually
    # aligned (as seen if one changes the semiglobal entry of question1c to
    # True)
    question1a = {"initial_gap_pen": 8, "end_gap_pen": 8, "semiglobal": False}
    question1b= {"initial_gap_pen": 4, "end_gap_pen": 4, "semiglobal": False}
    question1c = {"initial_gap_pen": 8, "end_gap_pen": 0, "semiglobal": False}

    print("Answer to question 1a:")
    question_1_print(seq1, seq2, 
                         question1a["initial_gap_pen"], 
                         end_gap_penalty = question1a["end_gap_pen"],
                         semiglobal = question1a["semiglobal"]
                         )

    print("Answer to question 1b:")
    question_1_print(seq1, seq2, 
                         question1b["initial_gap_pen"], 
                         end_gap_penalty = question1b["end_gap_pen"],
                         semiglobal = question1b["semiglobal"]
                         )

    print("Answer to question 1c:")
    question_1_print(seq1, seq2, 
                         question1c["initial_gap_pen"], 
                         end_gap_penalty = question1c["end_gap_pen"],
                         semiglobal = question1c["semiglobal"]
                         )
    print("\n\n")

    print("## ====================== QUESTION 3 ============================")
    
    # if you want to see the full details of each alignment as Q1, please
    # include an input of full_output = True when calling question_3()
    question_3(seq1, seq2)

    print("\n\n")

    print("## ====================== QUESTION 4 ============================")
    question_4(seq3, seq4, 5, 1)
    question_4(seq3, seq4, 5, 10)