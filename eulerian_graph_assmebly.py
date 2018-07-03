#!/usr/bin/env python

"""
Author: Stephen Coleman
"""
# Import statements
from difflib import SequenceMatcher
import warnings
from copy import deepcopy

def is_eulerian(graph):
    """Checks if a connected graph is Eulerian using the theorem:
    
            A connected graph is Eulerian if and only if each of its vertices
            is balanced

    graph: dict of lists, representing a connected graph
    """
    for key in graph.keys():
        out_degree = len(graph[key])

        in_degree = len([1 for paths in graph.values() if key in paths])

        if out_degree != in_degree:
            return False
    return True

def has_eulerian_path(graph):
    """Checks if a graph has an Eulerian path on the basis of the thoerem:

            A connected graph has an Eulerian path if and only if it contains 
            at most two semibalanced vertices and all other vertices are 
            balanced.

    graph: dict of lists, representing a connected graph
    """
    num_semiblanced = 0
    for key in graph.keys():
        out_degree = len(graph[key])

        in_degree = len([1 for paths in graph.values() if key in paths])

        if out_degree == in_degree:
            continue
        elif (out_degree - in_degree) ** 2 == 1 and num_semiblanced < 2:
            num_semiblanced += 1
        else:
            return False

    if num_semiblanced == 1:
        return False

    return True

def spectrum(l, s):
    """Returns a list of overlapping lmers in s

    l: int, length of lmers
    s: string, sequence of interest
    """
    if not isinstance(l, int):
        raise TypeError("l must be int")
    if l < 1:
        raise ValueError("l must be a positive non-zero integer")
    if not isinstance(s, str):
        raise TypeError("s must be a string")

    return [s[i : i + l] for i in range(0, len(s) - l + 1)]

def overlap(string1, string2):
    """Returns the maximum overlap between two strings

    string1: str
    string2: str
    """
    match = SequenceMatcher(None, string1, string2).\
              find_longest_match(0, len(string1), 0, len(string2))

    common_section = string1[match.a: match.a + match.size]
    return len(common_section)

def find_semi_balanced(graph):
    """Finds the semi balanced nodes in a graph

    graph: dict of lists, representing a connected graph
    """
    semibalanced_vertices = {}
    for key in graph.keys():
        out_degree = len(graph[key])
        
        in_degree = len([1 for paths in graph.values() if key in paths])

        if (out_degree - in_degree) ** 2 == 1:
            semibalanced_vertices[key] = [out_degree, in_degree, 
                                          out_degree - in_degree]
            
    return semibalanced_vertices

def find_longest_path(graph, start, end, path=[]):
    """Returns the longest path between two vertices on a graph

    graph: dict of lists, representing a connected graph
    start: key of graph, vertex in graph to begin at
    end: key of graph, vertex in graph to end at
    path: the path currently known between two vertices (for recursion 
    purposes, default = [])
    """
    path = path + [start]
    if start == end:
        return path
    if start not in graph:
        return None
    longest = None
    for node in graph[start]:
        if node not in path:
            newpath = find_longest_path(graph, node, end, path)

            if newpath:
                if not longest or len(newpath) > len(longest):
                    longest = newpath
    return longest

def find_eulerian_subcycles(graph):
    """ Returns a list of all eulerian cycles in a graph

    graph: dict of lists, representing a connected graph

    """
    if not is_eulerian(graph):
        raise Exception("graph is non-eulerian")

    start_point = list(graph.keys())[0]
    end_point = [key for key in graph.keys() if start_point in graph[key]][0]

    available_paths = deepcopy(graph)
    cycles = []

    while any(available_paths.values()):
        cycle = find_longest_path(available_paths, start_point, end_point)
        cycle += [start_point]
        cycles += [cycle]
        
        for i in range(len(cycle) - 1):
            available_paths[cycle[i]].remove(cycle[(i + 1) % len(cycle)])
        
        for vertex in cycle:
            if available_paths[vertex]:
                start_point = vertex
                end_point = [key for key in available_paths.keys() 
                             if start_point in available_paths[key]
                             ][0]
                break
    return cycles

def get_first_common_element(x,y):
    """ Returns first element from x that is common for both lists else None

    x: list
    y: list
    """
    for i in x:
        if i in y:
            return i

def merge_eulerian_cycles(list_eulerian_cycles):
    """Recursive function mergine list of eulerian cycles into one cycle

    list_eulerian_cycles: list of lists, each entry is an eulerian cycle. The 
    first entry should be connected to the second, second to the third, etc.
    """
    num_paths = len(list_eulerian_cycles)
    if num_paths == 1:
        return(list_eulerian_cycles[0])

    merge_point = get_first_common_element(list_eulerian_cycles[0],
                                           list_eulerian_cycles[1])
    if not merge_point:
        print("Not connected")
        return None

    merge_point_index_list_0 = list_eulerian_cycles[0].index(merge_point)
    merge_point_index_list_1 = list_eulerian_cycles[1].index(merge_point)

    merged_cycle = list_eulerian_cycles[0][ : merge_point_index_list_0] \
                   + list_eulerian_cycles[1][ : -1] \
                   + list_eulerian_cycles[0][merge_point_index_list_0 :]
    del list_eulerian_cycles[: min(2, len(list_eulerian_cycles))]
    list_eulerian_cycles = [merged_cycle] + list_eulerian_cycles
    return merge_eulerian_cycles(list_eulerian_cycles)

def check_all_list_entries_equal_length(my_list):
    """Raises an error if the entries in a list are not all of common length
    
    my_list: list
    """
    it = iter(my_list)
    the_len = len(next(it))
    if not all(len(l) == the_len for l in it):
        raise ValueError('not all entries have the same length.')

def graph_from_sequences(sequences):
    """Returns a graph of (l-1)-mers from a sequence of lmers

    sequences: list of strings, each string must be an lmer of common l
    """
    check_all_list_entries_equal_length(sequences)
    l = len(sequences[0]) # should check that all sequences are the same length
    vertices = {seq : spectrum(l - 1, seq) for seq in sequences}
    graph = {}

    # Use a set to ensure only unique (l-1) mers are recorded
    vertices_values = {vertex for sublist in vertices.values() 
                       for vertex in sublist
                       }

    for vertex in vertices_values:
        graph[vertex] = [edge[1:] for edge in vertices.keys()
                         if vertex == vertices[edge][0]
                         ]

    return graph


def find_eulerian_cycle(graph):
    """Returns the Eulerian cycle (if present) for a connected graph

    graph: dict of lists, representing a connected graph
    """
    subcycles = find_eulerian_subcycles(graph)
    return merge_eulerian_cycles(subcycles)

def find_eulerian_path(graph):
    """Returns the eulerian path through a graph as an ordered list of vertices

    graph: dict of lists, representing a connected graph
    """
    
    if is_eulerian(graph):
        return find_eulerian_cycle(graph) # do not need the end point
    if not has_eulerian_path(graph):
        warnings.warn("no Eulerian path present in graph.")
        return None

    local_graph = deepcopy(graph)

    # Find the semibalanced vetices - these will be our start/end points
    semibalanced_vertices = find_semi_balanced(graph)

    # Order them based on which has a lower relative out degree to in degree
    semibalanced_order = sorted(semibalanced_vertices.items(), 
                                key=lambda x: x[1][2]
                                )

    # Add an edge from the end point (that with lower out degree than in 
    # degree) to the start point
    local_graph[semibalanced_order[0][0]] += [semibalanced_order[1][0]]

    # Now that the vertices are all balanced, treat as an eulerian cycle
    eulerian_cycle = find_eulerian_cycle(local_graph)

    # The actual start point of our path and its index within the pseudo-cycle
    start_point = semibalanced_order[1][0]
    start_index = eulerian_cycle.index(start_point)

    # Reorder the cycle and drop the original cycle's end entry (as is a 
    # duplicate of the start)
    eulerian_path = eulerian_cycle[start_index : ] \
                     + eulerian_cycle[1 : start_index]

    return eulerian_path

def print_graph(graph):
    """Prints graph based on the request from the assignment doc

    graph: dict of lists, representing a graph    
    """
    for k, v in graph.items():
        print(k, v)


def sequence_from_eulerian_path(eulerian_path):
    """Returns a sequence from an eulerian path of lmers with overlap of (l-1)

    eulerian_path: output of find_eulerian_path for a graph of sequences
    """

    sequence = eulerian_path[0]
    if len(eulerian_path) == 1:
        return sequence

    for lmer in eulerian_path[1 :]:
        sequence += lmer[-1]
    return sequence

def printed_output(graph, times_check_the_same = 3, combine_sequence = False):
    """Printed output for questions 1 - 4 of assignment

    graph: dict of lists, representing a graph    
    """
    print("\nGraph:")
    print_graph(graph)
    print("\nIs Eulerian: {}".format(is_eulerian(graph)))
    print("Contains an Eulerian path: {}". \
          format(has_eulerian_path(graph))
          )
    print("Eulerian path:\n{}".format(find_eulerian_path(graph)))

    paths = set()
    always_the_same = True
    for i in range(times_check_the_same):
        new_path = find_eulerian_path(graph)
        paths.add(tuple(new_path))
    if len(paths) != 1:
        always_the_same = False
    print("Path always the same: {}\n".format(always_the_same))

    if combine_sequence:
        print("Sequence: {}\n". \
              format(sequence_from_eulerian_path(find_eulerian_path(graph)))
              )

if __name__ == "__main__":

    # GRAPH FROM FIG 8.22
    graph_822 = {'A':['B'],'B':['C'],'I':['H'],'H':['F'],'F':['G','E'],\
        'C':['I','J'],'G':['A'],'E':['J'],'J':['F','D'],'D':['C']}



    # A SLIGHTLY BIGGER GRAPH, NEEDED FOR Q8

    bigger_graph = {1:[2], 2:[3], 3:[7,1],\
        4:[5,10],5:[6],6:[7],7:[8,9,11],\
        8:[4],9:[3,4],\
        10:[9],11:[12],12:[7]}

    # SPECTRUM FROM FIG 8.20
    s = ['ATG','TGG','TGC','GTG','GGC','GCA','GCG','CGT']

    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file) 
    
    print("=== Fig. 8.22 ===================================================")
    printed_output(graph_822)

    print("=== Sequence ====================================================")
    seq_graph = graph_from_sequences(s)
    printed_output(seq_graph, combine_sequence = True)
    print("Eulerian path present as there are exactly two semi-balanced \
           \nvertices present, and all others are balanced and connected.\n")

    print("=== Bigger graph ================================================")
    printed_output(bigger_graph)