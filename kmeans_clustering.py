#!/usr/bin/env python
"""
Author: Stephen Coleman
Implementation of the k-means clustering algorithm
"""

# import statements
import random, sys, warnings
from math import sqrt, isclose
from numpy import argmin, mean, sum, subtract, std, exp
from numpy.random import uniform, rand
import matplotlib.pyplot as plt
from scipy.stats import mode
from copy import deepcopy

def csv_parser(lines):
    """Return list of point coordinates as [[x1,y1,z1,...],[x2,y2,z2,...]]
    
    lines: open file or list of lines. Expected format:
        The file has a single header line. 
        Each line contains the coordinates for one data point, starting
        with a label. A data point can be specified in arbitrary dimensions.

    Output: List of lists with the coordinates of the points.
    This function does not capture the labels of the data points. In case
    they are needed later, this function should be adjusted. 
    """ 

    data_points = []
    for line in lines:
        items = line.strip().split(",")
        try: #will fail on header line in file
            data_points.append(list(map(float, items[1:]))) #first item is the label
        except ValueError: #must be the header
            continue
    return data_points

def euclidean(point1, point2):
    """Returns the euclidean distance between vectors x and y

    point1: list of spatial coordinates
    point2: list of spatial coordinates of equal dimension to point1
    """
    if len(point1) != len(point2):
        raise ValueError("point1 must have the same dimension as point2")
    
    dist = sqrt(sum(subtract(point1, point2) ** 2))

    return dist

def initialise_centres(points, k):
    """Returns k random points from points to use as initial centres

    data: list of spatial coordinates
    k: number of centres to use
    """
    if len(points) < k:
        raise ValueError("k exceeds the number of points present." \
                          "\nPlease use a value less than or equal to {}"\
                         .format(len(points))
                         )
    return random.sample(points, k)

def initialise_random_centre(data, k):
    """
    """
    infimum = find_minimum(data)
    sup = find_maximum(data)
    dim = len(sup)
    centres_unit = rand(k, dim)
    coefficients = uniform(low = min(infimum), high = max(sup), size = k)
    centres = []
    for i, mult in enumerate(coefficients):
        centres += [centres_unit[i] * mult]
    # print(centres)
    return centres

def find_minimum(data):
    """Returns the lowest point of the hyperrectangle containing all data

    data: list of spatial coordinates
    """
    d = len(data[0])
    vector_components = list(zip(*data))

    minimum = [0] * d

    for i, dimension in enumerate(vector_components):
        minimum[i] = min(dimension)

    return(minimum)

def find_maximum(data):
    """Returns the highest point of the hyperrectangle containing all data

    data: list of spatial coordinates
    """
    d = len(data[0])
    vector_components = list(zip(*data))

    maximum = [0] * d

    for i, dimension in enumerate(vector_components):
        maximum[i] = max(dimension)

    return(maximum)

def order_data(data):
    """Returns data sorted on order in Euclidean space

    data: list of spatial coordinates
    """
    minimum = find_minimum(data)
    distances = []
    for point in data:
        distances += [euclidean(point, minimum)]

    ordered_data = sorted(list(zip(distances, data)))
    ordered_data = list(zip(*ordered_data))[1]
    return list(ordered_data)


def split_along_median(ordered_data):
    """Returns the median of ordered data

    ordered_data: ordered list / tuple of coordinates
    """
    split_index = (len(ordered_data) // 2) + 1
    return [ordered_data[: split_index], ordered_data[split_index :]]

def kd_tree(data, max_points = 1):
    """Returns a list of nodes clustered based on Euclidean distance

    data: list of spatial coordinates
    max_points: int; upper bound on number of points to include in the final
    node (default = 1)
    """
    ordered_data = order_data(data)
    tree = [ordered_data]

    i = 0
    while any(len(node) > max_points for node in tree):
        node = tree[i]
        while len(node) > max_points:
            tree = tree[:i] + split_along_median(node) + tree[i + 1:]
            node = split_along_median(node)[0]
        i += 1
    return tree

def nearest_centre(point, centres):
    """Returns the nearest centre for a given point from a set of centres

    point: list, spatial coordinates
    centres: list of spatial coordinates
    """
    distances = []
    for centre in centres:
        distances += [euclidean(point, centre)]

    # return centres[argmin(distances)]
    return argmin(distances)

def assign_classes(data, centres):
    """Returns a list classifying the points in data based on centres

    data: list of spatial coordinates
    centres: list of cluster centres
    """
    groups = []
    for point in data:
        groups += [nearest_centre(point, centres)]
    return groups

def update_centres_old(data, groups):
    """Returns the centres of the groups associated with data

    data: list of spatial coordinates
    groups: list of classificaitons, where ith entry corresponds with the 
    ith entry of data
    """
    centres = []
    drop = []
    k = set(groups)
    grouped_data = list(zip(groups, data))
    for i in k:
        cluster = [point[1] for point in grouped_data if point[0] == i]
        centres += [tuple(mean(cluster, axis = 0))]
    
    return centres

def update_centres(data, groups, old_centres):
    """Returns the centres of the groups associated with data

    data: list of spatial coordinates
    groups: list of classificaitons, where ith entry corresponds with the 
    ith entry of data
    """
    centres = []
    present = []
    k = set(groups)
    grouped_data = list(zip(groups, data))
    for i in k:
        cluster = [point[1] for point in grouped_data if point[0] == i]
        centres += [tuple(mean(cluster, axis = 0))]
        present += [old_centres[i]]

    if len(present) != len(old_centres):
        warnings.warn("No points assigned to some clusters: dropping clusters\
                       \nUpdated number of clusters = {}.".format(len(present))
                      )
    return centres, present

def kmeans(data, k, 
           threshold = 0.00001,
           verbose = False, 
           restrict_initial_values = True
           ):
    """Returns the cluster centres and classification list for data

    data: list of spatial coordinates
    k: int; number of clusters to use
    threshold: threshold for movement of clusters between iteraitons 
    (default = 0.00001)
    verbose: bool; instructs printing information upon each iteration 
    (default = False)
    restrict_initial_values: bool; instructs kmeans to limit initial centres
    to random points within data or random points within the same range as 
    data (default = True)
    """
    if restrict_initial_values:
        centres = initialise_centres(data, k)
    else:
        centres = initialise_random_centre(data, k)
    
    movement = float("Inf")
    i = 0
    while movement > threshold:
        groups = assign_classes(data, centres)

        if i == 0:
            plot_kmeans(data, centres, groups, marker = "x")
        new_centres, centres = update_centres(data, groups, centres)
 
        # Manhattan metric
        # movement = sum(abs(subtract(new_centres, centres)))

        # Euclidean metric
        movement = sum([euclidean(new, centres[i]) 
                        for i, new in enumerate(new_centres)
                        ]
                       )

        centres = new_centres
        if verbose:
            print("Iteration {}: movement = {:.2f}".\
                  format(i, movement)
                  )
            for j, centre in enumerate(centres):
                coords_string = ", ".join("{:.2f}".format(point) 
                                          for point in centre
                                          )

                print("Centre {}: ({})".format(j, coords_string))
            print("\n")
        i += 1

    # print("Convergence in {} iterations.".format(i))

    return centres, groups, i

def plot_kmeans(data, centres, groups, marker = "+"):
    """Plots kmeans data in first two dimensions

    data: list of spatial coordinates
    centres: list of cluster centres
    groups: list of assigned classes for points in data
    """
    plot_data = list(zip(*data))
    plot_centres = list(zip(*centres))
    k = len(centres)

    with plt.style.context("ggplot"):
        plt.scatter(plot_data[0], plot_data[1], c = groups)
        plt.scatter(plot_centres[0], plot_centres[1], 
                    c = range(k),
                    marker = marker
                    )

def WGCSS_c(cluster_data, centre):
    """Returns the within-cluster sum of squares for a single cluster

    cluster_data: list of points associated with a single cluster
    centre: the coordinates of the associated centre
    """
    count = len(cluster_data)

    distances = [euclidean(point, centre) ** 2 for point in cluster_data]
    return sum(distances) / count

def WGCSS(data, centres, groups):
    """Returns the total within-cluster sum of squares for all clusters

    data: list of spatial coordinates
    centres: list of cluster centres
    groups: list of assigned classes for points in data
    """
    cluster_WGCSS_c = [WGCSS_c(find_cluster_data(data, groups, i), centre) 
                       for i, centre in enumerate(centres)
                       ]

    return sum(cluster_WGCSS_c)

def find_cluster_data(data, groups, group_number):
    """Returns the subset of data assigned to cluster group_number

    data: list of spatial coordinates
    groups: list of assigned classes for points in data
    group_number: int; class of interest
    """
    return [point for i, point in enumerate(data) if groups[i] == group_number]

def BGCSS(data, centres, groups):
    """Returns the between cluster sum of squares
    
    data: list of spatial coordinates
    centres: list of cluster centres
    groups: list of assigned classes for points in data
    """
    data_centre = mean(data, axis = 0)
    counts = []
    for i in range(len(centres)):
        counts += [len(find_cluster_data(data, groups, i))]

    distances = [euclidean(centre, data_centre) ** 2 for centre in centres]

    return sum([counts[i] * distances[i] for i in range(len(distances))])

def criterion(data, centres, groups):
    """Returns the ratio of WGCSS to BGCSS for clusterings within data

    data: list of spatial coordinates
    centres: list of cluster centres
    groups: list of assigned classes for points in data
    """
    num = WGCSS(data, centres, groups)
    denom = BGCSS(data, centres, groups)
    return  num / denom


def print_output(filename, k, 
                 verbose = False, 
                 restrict_initial_values = True,
                 num_coord_dim_to_print = 2,
                 variability_incl = False,
                 num_iter = 101,
                 plot = True,
                 print_groups = False,
                 ):
    """Prints the requested output of kmeans clustering
    
    filename: str; name of csv file containing data
    verbose: bool; instructs printing of progress at each iteraiton of 
    clustering (default is False)
    restrict_initial_values: bool; restricts initial centres to points within
    data contained in filename (default is True)
    num_coord_dim_to_print: int; number of dimensions of coordinates to print
    if larger than dimensionality of point, defaults to dimensionality of 
    point (defualt is 2)
    variability_incl: bool; instructs recording of variability of W
    num_iter: int; only relevant if variability_incl is True, number of 
    iterations to record W over
    plot: bool; instructs plotting of clusters
    print_groups: bool; instructs printing of classifications
    """

    best_W = float("Inf")
    best_k = float("Inf")

    if isinstance(k, int):
        k = [k]

    with open(filename) as f:
        data = csv_parser(f)

    print("\nDataset: {}".format(filename))

    for num_clusters in k:
        
        print("k: {}".format(num_clusters))

        if variability_incl:
            centres, mean_W, sd_W, min_W, num_iterations, groups, W_range = \
                lots_kmeans(data, num_clusters, 
                            verbose = verbose,
                            restrict_initial_values = restrict_initial_values,
                            num_iter = num_iter
                            )

            print("Across {} iterations:".format(num_iter))
            # print("W: {:.3f} +/- {:.3f}".format(mean_W, 2 * sd_W))
            print("Mean W: {:.3f}".format(mean_W))
            print("Standard deviation W: {:.3f}".format(sd_W))
            print("Min W: {:.3f}".format(min_W))

            print("For mode cluster:")
            if min_W < best_W:
                best_W = min_W
                best_k = num_clusters
            
        else:
            centres, groups, num_iterations = kmeans(data, num_clusters, 
                                                     verbose = verbose,
                                                     restrict_initial_values =\
                                                     restrict_initial_values
                                                     )

            
            W = criterion(data, centres, groups)

            print("W: {:.3f}".format(W))

            if W < best_W:
                best_W = W
                best_k = num_clusters

        print("Convergence in {} iterations".format(num_iterations))


        if print_groups:
            print(groups)

        for i, centre in enumerate(centres):
            # print(centre)
            centre_string = coordinate_string(centre, 
                                              num_dim = num_coord_dim_to_print
                                              )

            print("Centre group {}:\t{}".\
                    format(i, centre_string))

            num_members = len([member for member in groups if member == i])
            print("Number members: {}".format(num_members))

        print("\n")
        
        if plot:
            with plt.style.context("ggplot"):
                plot_kmeans(data, centres, groups)
                plt.suptitle("K means clustering")
                plt.title("\nData: {}, K = {}".format(filename, num_clusters), 
                           fontsize=9, 
                           color='grey',
                           style='italic')

            plt.show()

    print("\nAcross {}, optimal k = {} with W = {:.3f}".\
                format(filename, best_k, best_W)
              )

def lots_kmeans(data, k,
                verbose = False,
                restrict_initial_values = True,
                num_iter = 101
                ):

    """Returns the mode centres and variability of W for kmeans clustering

    data: list of spatial coordinates
    k: int; number of clusters
    verbose: bool; instructs printing of progress at each iteraiton of 
    clustering (default is False)
    restrict_initial_values: bool; restricts initial centres to points within
    data contained in filename (default is True)
    num_iter: number of iterations for carrying out kmeans
    """
    W_range = []
    centres_range = []
    num_iter_range = []
    groups_range = []
    for i in range(num_iter):
        centres, groups, num_iterations = kmeans(data, k, 
                                                 verbose = verbose,
                                                 restrict_initial_values = \
                                                     restrict_initial_values
                                                 )
        W_range += [criterion(data, centres, groups)]
        centres_range += [centres]
        num_iter_range += [num_iterations]
        groups_range += [groups]

    # Below we record the mode of the centres across the iterations
    # We include the try loop to avoid multiple modes
    mode_centre = None
    while mode_centre is None:
        try:
            # most frequently recurring centres
            mode_centre = mode(centres_range)

        except ValueError as Err:
            warnings.warn("Multiple modes, removing an option")
            centres_range = centres_range[:-1]

    mean_W = mean(W_range)
    sd_W = std(W_range)
    min_W = min(W_range)

    mode_centre = mode_centre[0][0]
    
    # Find the grouping and number of iterations from the mode centre
    for i, centre in enumerate(centres_range):
        same = []
        for j, point in enumerate(centre):
            for k, dim in enumerate(point):
                same += [isclose(dim, mode_centre[j][k])]
                if all(same):
                    final_grouping = groups_range[i]
                    final_num_iter = num_iter_range[i]
                    break
    
    return mode_centre, mean_W, sd_W, min_W, final_num_iter, final_grouping, W_range   

def coordinate_string(point, num_dim = None, precision = 2):
    """Converts a point to a string of given precision and dimensions

    point: lsit / tuple of spatial coordinates (as floats)
    num_dim: number of dimensions to print
    precision: precision of each point
    """
    if num_dim > len(point):
        warnings.warn("Dimensions requested exceed dimensionality of point\
                       \nDefaulting to point dimensionality")
        num_dim = len(point)
    if num_dim is None:
        num_dim = len(point)
    coords_string = "("
    for dim in point[:num_dim]:
        curr_coord = "{:.{precision}f},".format(dim, precision = precision)
        coords_string += curr_coord
    coords_string = coords_string[:-1]
    coords_string += ")"
    return coords_string


def bootstrap(data):
    """Create new sample via bootstrapping

    data: list of spatial coordinates (as ints / floats)
    """
    new_data = [data[random.randint(0, len(data) - 1)] for point in data]
    return new_data

def subsample(data, subsample_size = None):
    """Returns a subsample of data without replacement

    data: list of spatial coordinates (as ints / floats)
    subsample_size: int; if left as None uses (1-e**-1 ~ 0.63) * len(data)
    (default is None)
    """
    if subsample_size is None:
        subsample_size = int((1 - exp(-1)) * len(data))

    if subsample_size == len(data):
        return data

    loc_data = deepcopy(data)

    new_data = []
    for i in range(subsample_size):
        index_new_point = random.randint(0, len(data) - 1)
        new_point = loc_data[index_new_point]

        new_data += [new_point]
        
    return new_data

def W_over_resample(filename, k, method,
                    verbose = False, 
                    restrict_initial_values = True,
                    num_coord_dim_to_print = 2,
                    num_iter = 11,
                    resampling_iter = 10
                    ):
    """Prints the requested output of kmeans clustering
    
    filename: str; name of csv file containing data
    k: list of ints; number of clusters
    method: string; one of bootstrap or subsample
    verbose: bool; instructs printing of progress at each iteraiton of 
    clustering (default is False)
    restrict_initial_values: bool; restricts initial centres to points within
    data contained in filename (default is True)
    num_coord_dim_to_print: int; number of dimensions of coordinates to print
    if larger than dimensionality of point, defaults to dimensionality of 
    point (defualt is 2)
    num_iter: int; only relevant if variability_incl is True, number of 
    iterations to record W over
    resampling_iter: int; number of times to resample data
    """

    if method not in ["bootstrap","subsample"]:
        raise ValueError("method must be one of bootstrap or subsample")

    W_tot_range = []

    with open(filename) as f:
        main_data = csv_parser(f)

    for resamp_iter in range(resampling_iter):
        
        if method == "bootstrap":

            data = bootstrap(main_data)
        else:
            data =subsample(main_data)

        centres, mean_W, sd_W, min_W, num_iterations, groups, W_range = \
            lots_kmeans(data, k, 
                        verbose = verbose,
                        restrict_initial_values = \
                            restrict_initial_values,
                        num_iter = num_iter
                        )

        W_tot_range += [W_range]

    return W_range

def question8(k, method, filename, 
              verbose = False, 
              restrict_initial_values = True,
              num_coord_dim_to_print = 2, 
              resampling_iter = 10):
    """Output for question 8


    """
    for num_clusters in k:
        W_range = W_over_resample(filename, num_clusters, method,
                                  verbose = verbose,
                                  restrict_initial_values = \
                                    restrict_initial_values,
                                  num_coord_dim_to_print = 24,
                                  resampling_iter = 10
                                  )
        
        min_W = min(W_range)
        mean_W = mean(W_range)
        sd_W = std(W_range)

        print("k: {}".format(num_clusters))
        print("Mean W: {:.3f}".format(mean_W))
        print("Standard deviation W: {:.3f}".format(sd_W))
        print("Min W: {:.3f}".format(min_W))

if __name__ == "__main__":

    # the code below should produce the results necessary to answer
    # the questions. In other words, if we run your code, we should see 
    # the data that you used to answer the questions.
      
    verbose = False
    restrict_initial_values = True

    print("=== Question 1: 2dtest.csv ======================================")
    filename = "2dtest.csv"
    k = 3
    print_output(filename, k, 
                 verbose = verbose,
                 restrict_initial_values = restrict_initial_values
                 )

    print("=== Question 2: LargeSet_1.csv ==================================")
    filename = "LargeSet_1.csv"
    k = range(2, 7)
    for i in range(5):
        print_output(filename, k, 
                     verbose = verbose,
                     restrict_initial_values = restrict_initial_values,
                     plot = False
                     )

    print("=== Question 3: LargeSet_2.csv ==================================")
    filename = "LargeSet_2.csv"
    k = [2] * 5
    print_output(filename, k, 
                 verbose = verbose,
                 restrict_initial_values = restrict_initial_values
                 )

    print("=== Question 4: LargeSet_1.csv  variability =====================")
    filename = "LargeSet_1.csv"
    k = range(2, 7)
    print("\n")
    print_output(filename, k, 
                 verbose = verbose,
                 restrict_initial_values = restrict_initial_values,
                 variability_incl = True,
                 num_coord_dim_to_print = 24,
                 plot = False,
                 print_groups = True
                 )

    print("=== Question 8: Resampling ======================================")
    filename = "LargeSet_1.csv"
    k = range(2, 7)
    method = "bootstrap"
    print("\n")
    print("Sampling method: {}".format(method))
    question8(k, method, filename)
    
    method = "subsample"
    print("\n")
    print("Sampling method: {}".format(method))
    question8(k, method, filename)