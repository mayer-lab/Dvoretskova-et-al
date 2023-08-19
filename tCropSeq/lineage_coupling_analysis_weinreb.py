import codecs
from contextlib import closing
import requests
import csv
import itertools
import statistics
import random
import collections
import copy
import math
import pprint as pp
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import argparse


def parse_args():
    """Parse the arguments specified by the user."""
    vmin = 0.0
    vmax = 3.0
    num_iterations = 100
    csv_data_url = ('https://keeper.mpdl.mpg.de/lib'
                    '/011bef7d-af64-4883-9135-afebe1e25e1e/file/results/TIJ'
                    '/results/lb_pool/TIJ_HARMONY_Clusters.csv?dl=1')
    csv_data_file = 'Datasets/TIJ_HARMONY_Clusters.csv'
    csv_metric_values_matrix_output_file = (''
            'Results/metric_values_per_state_pair_matrix.csv')
    csv_lineage_coupling_scores_matrix_output_file = (''
                                'Results/lineage_coupling_scores_matrix.csv')
    parser = argparse.ArgumentParser(description="Perform"
                    " the Lineage Analysis derived from Weinreb et al. (2020)"
                    " for a given input.")
    parser.add_argument("-N", "--num_iterations", metavar='N'
                        , type=int, nargs=1
                        , help="The number of iterations"
                                " (of clones sampling) to be made"
                                " (default {}).".format(num_iterations))
    parser.add_argument("-l", "--csv_url", metavar='csv_data_url'
                , type=str, nargs=1
                , help="The csv url containing the data"
                        " (default '{}').".format(csv_data_url))
    parser.add_argument("-f", "--csv_file", metavar='csv_data_file'
                        , type=str, nargs=1
                        , help="The csv file containing the data"
                                " (default '{}').".format(csv_data_file))
    parser.add_argument("-M"
                        , "--csv_metric_values_matrix_output_file"
                        , metavar='csv_metric_values_matrix_output_file'
                        , type=str, nargs=1
                        , help="The output csv file where the "
                                "metric values per state pair will be written"
                                " (note that the folder must exist)"
                                " (default '{}').".format(
                                    csv_metric_values_matrix_output_file)
            )
    parser.add_argument("-L"
                , "--csv_lineage_coupling_scores_matrix_output_file"
                , metavar='csv_lineage_coupling_scores_matrix_output_file'
                , type=str, nargs=1
                , help="The output csv file where the lineage coupling"
                    "scores will be written (note that the folder must exist)"
                    " (default '{}').".format(
                            csv_lineage_coupling_scores_matrix_output_file)
            )
    parser.add_argument("-u", "--clustermap_scale_min", metavar='scale_min'
                        , type=float, nargs=1
                        , help="The minimum number of the scale used"
                                " to plot the clustermap (default {})."
                                                                .format(vmin))
    parser.add_argument("-v", "--clustermap_scale_max", metavar='scale_max'
                        , type=float, nargs=1
                        , help="The maximum number of the scale used"
                                " to plot the clustermap (default {})."
                                                                .format(vmax))
    args = parser.parse_args()
    if(args.num_iterations is not None):
        num_iterations = args.num_iterations[0]
    # The url has priority over the path
    if(args.csv_url is not None):
        csv_data_url = args.csv_url[0]
        csv_data_file = None
    elif(args.csv_file is not None):
        csv_data_url = None
        csv_data_file = args.csv_file[0]
    if(args.csv_lineage_coupling_scores_matrix_output_file is not None):
        csv_lineage_coupling_scores_matrix_output_file = (
                        args.csv_lineage_coupling_scores_matrix_output_file[0]
                    )
    if(args.csv_metric_values_matrix_output_file is not None):
        csv_metric_values_matrix_output_file = (
                    args.csv_metric_values_matrix_output_file[0]
                )
    if(args.clustermap_scale_min is not None):
        vmin = args.clustermap_scale_min[0]
    if(args.clustermap_scale_max is not None):
        vmax = args.clustermap_scale_max[0]
    return (num_iterations
            , csv_data_url, csv_data_file
            , csv_lineage_coupling_scores_matrix_output_file
            , csv_metric_values_matrix_output_file, vmin, vmax)

def extract_rows(csv_reader, csv_data_clone_cells):
    """
    Extract relevant row values from a csv.DictReader() object.

    Obs:
        In some datasets there may be a column, 'ident_name'
        which stores the names of the cell states.
        That column is used for parts that should output something
        about the cell states (e.g. metric values in the csv).
        However, in all datasets used for the analyses in the paper,
        there wasn't such column, so the column 'ident' itself is used
        as the column name.
    """
    row_tmp = {}
    for i, row in enumerate(csv_reader):
        # Store the register in a new dedicated list.
        # Discard all registers that had an empty value
        # for any of the fields
        # Also, for some reason, if cloneID is converted into a string
        # rather than a integer, then other, unrelated
        # floating computations become a little bit unstable 
        # (e.g. the metric values begin to differ in
        # their latter decimals when computed with
        # the exact same parameters)
        row_tmp["cloneID"] = int(row["cloneID"])
        if(row_tmp["cloneID"] == ""):
            print("Register number {}.".format(i))
            print("cloneID is empty."
                    .format(row["cloneID"])
                    )
            continue
        row_tmp["ident_name"] = str(row["ident"])
        if(row_tmp["ident_name"] == ""):
            print("Register number {}.".format(i))
            print("ident_name is empty.\n"
                    .format(row["ident"])
                    )
            continue
        row_tmp["ident"] = str(row["ident"])
        if(row_tmp["ident"] == ""):
            print("Register number {}.".format(i))
            print("ident is empty.\n"
                .format(row["ident"])
                )
            continue
        csv_data_clone_cells.append(row_tmp.copy())
    return csv_data_clone_cells

def parse_csv_data(csv_data_url, csv_data_file):
    """
    Return the relevant parsed data from a selected source.

    The relevant data are constituted by:
        1) Data of cells that belong to a clone,
        without including the name of the cluster.
        2) The names of each cluster.
    """
    csv_data_clone_cells = []
    # Before calling this function, if the selected source was
    # a csv file, then csv_data_url was set to None. Otherwise,
    # the url source has priority.
    if(csv_data_url is not None):
        response = requests.get(csv_data_url)
        with closing(requests.get(csv_data_url, stream=True)) as req:
            csv_reader = csv.DictReader(codecs.iterdecode(
                                                req.iter_lines(), 'utf-8'))
            csv_data_clone_cells = extract_rows(csv_reader
                                                    , csv_data_clone_cells)
    else:
        with open(csv_data_file, newline='') as csvfile:
            csv_reader = csv.DictReader(csvfile)
            csv_data_clone_cells = extract_rows(csv_reader
                                                    , csv_data_clone_cells)
            csvfile.close()
    return csv_data_clone_cells

def compute_clones_distribution(csv_data_clone_cells):
    """Return a counter of cells per clone."""
    clones = collections.Counter()
    for row in csv_data_clone_cells:
        cloneID = row["cloneID"]
        clones[cloneID] += 1
    return clones

def compute_clusters_distribution(csv_data_clone_cells):
    """Return a counter of cells per cluster, and their names."""
    clusters = collections.Counter()
    clusters_names = {}
    for row in csv_data_clone_cells:
        clusterID = row["ident"]
        cluster_name = row["ident_name"]
        clusters[clusterID] += 1
        clusters_names[clusterID] = cluster_name
    return clusters, clusters_names

def generate_cluster_pairs(clusters_list):
    """
    Return a dictionary-valued dictionary indexed by the cluster pairs.

    The indices are 2-tuples, whose elements are the clusterIDs.
    Despite this, each pair is conceptually considered to be
    an unordered pair, and hereby appearing only once.
    This dictionary will be expanded to be used for other
    purposes as well.
    """
    cluster_pairs = {}
    for cluster_pair in itertools.combinations_with_replacement(
                                                    clusters_list, 2):
        cluster_pairs[cluster_pair] = {}
    return cluster_pairs

def compute_num_cells_per_clones_per_cluster(csv_data_clone_cells):
    """
    Return a multi-level counter with different cell counts per level.

    Each count is stored in a field named "total".
    What is counted on each level:
        1) Total number of cells.
        2) Total number of cells per cluster.
        3) Total number of cells per clone per cluster.
    """
    num_cells = collections.Counter()
    num_cells["total"] = 0
    for cell in csv_data_clone_cells:
        cluster = cell["ident"]
        clone = cell["cloneID"]
        num_cells["total"] += 1
        if(cluster not in num_cells):
            num_cells[cluster] = collections.Counter()
            num_cells[cluster]["total"] = 0
        num_cells[cluster]["total"] += 1
        if(clone not in num_cells[cluster]):
            num_cells[cluster][clone] = collections.Counter()
            num_cells[cluster][clone]["total"] = 0
        num_cells[cluster][clone]["total"] += 1
    return num_cells

def compute_total_num_cells_of_shared_clone_and_cluster_pair(
    clone, cluster_pair, num_cells):
    """
    Return the number of cells in a cluster pair and its shared clone.

    For a clone and a cluster pair,
    if the clone was "shared" by the cluster pair (according to the
    corresponding definition in the Methods section of the paper),
    count the cells that belonged to the clone
    and also were assigned to one of the two clusters.
    If the clone was not shared between the clusters, return 0.
    """
    total_num_cells_of_shared_clone_and_cluster_pair = 0
    # The number of cells of the shared clone of each cluster is added
    for cluster in cluster_pair:
        # If no cells of a clone were assigned to any of the clusters
        if(clone not in num_cells[cluster]
                    or num_cells[cluster][clone]["total"] < 1):
            # Then that clone isn't "shared" by the pair
            return 0
        # If the cluster pair is composed by the same cluster
        if(cluster_pair[0] == cluster_pair[1]):
            # Then the total number of cells of the shared clone
            # is just that of that cluster
            total_num_cells_of_shared_clone_and_cluster_pair = (
                                        num_cells[cluster][clone]["total"])
            break
        else:
            total_num_cells_of_shared_clone_and_cluster_pair += (
                                        num_cells[cluster][clone]["total"])
    # If there were less than 2 cells of a clone
    # that were assigned to any of the clusters, then
    # that clone isn't "shared" between them
    if(total_num_cells_of_shared_clone_and_cluster_pair < 2):
        return 0
    else:
        return total_num_cells_of_shared_clone_and_cluster_pair

def compute_metric_value_single_pair(num_cells, cluster_pair
    , clones):
    """
    Return the metric value of a cluster pair.

    (According to the definition
    given in the Methods section of the paper).
    """
    metric_value = 0
    # Check the case where
    # a cluster doesn't have any cells assigned to it,
    # just because this method is called when doing the samples,
    # and within one of them, this case may arise.
    # Using the observed data this shouldn't happen, though,
    # given that cluster_pairs in the calling method
    # is generated so as abscent clusters are omitted.
    if(num_cells[cluster_pair[0]] == 0 or num_cells[cluster_pair[1]] == 0):
        return 0
    # For each clone that had cells assigned to any of the states
    for clone in set(itertools.chain(num_cells[cluster_pair[0]].keys()
                                    , num_cells[cluster_pair[1]].keys())
                    ):
        # Note that in the previous operation, as num_cells stores
        # the total value in a separated field, in the same level that
        # the clones, the "total" field is also included
        if(clone == "total"):
            # So skip it
            continue
        # Only the cells that belonged to shared clones are added
        total_num_cells_of_shared_clone_and_cluster_pair = (
                compute_total_num_cells_of_shared_clone_and_cluster_pair(
                                            clone, cluster_pair, num_cells))
        # Compute and store the relative value with respect to the
        # total number of cells in the (shared) clone
        metric_value += (total_num_cells_of_shared_clone_and_cluster_pair
                                                            / clones[clone])
    return metric_value

def compute_observed_metric_value_all_pairs(num_cells, cluster_pairs
    , clones):
    """
    Return the metric value of each cluster pair.

    It is stored in a dictionary field, indexed by the cluster pair.
    """
    for cluster_pair in cluster_pairs:
        cluster_pairs[cluster_pair]["metric_value"] = (
                    compute_metric_value_single_pair(
                                num_cells, cluster_pair, clones))
    return cluster_pairs

def compute_sum_metric_value_cluster(cluster, cluster_pairs):
    """
    Return the sum of metric values of a cluster.

    This sum is the one obtained
    when the metric values of all cluster pairs
    in which the cluster was present are added together.
    It would be the sum of a given row (or column) of the Oij matrix.
    """
    sum_metric_value_by_cluster = 0
    for cluster_pair in cluster_pairs:
        if(cluster in cluster_pair):
            sum_metric_value_by_cluster += (
                        cluster_pairs[cluster_pair]["metric_value"])
    return sum_metric_value_by_cluster

def compute_sum_metric_value_all_clusters(cluster_pairs):
    """
    Return the total sum of all metric values.

    This sum is the one obtained
    when the metric values of all cluster pairs are added together.
    It would be the sum of the lower triangular Oij matrix
    and the diagonal.
    """
    sum_metric_value = 0
    for cluster_pair in cluster_pairs:
        sum_metric_value += cluster_pairs[
                                        cluster_pair]["metric_value"]
    return sum_metric_value

def compute_sum_metric_value_same_clusters(cluster_pairs):
    """
    Return the total sum of all metric values.

    This sum is the one obtained when the metric values of all
    cluster pairs composed by the same cluster are added together.
    It would be the sum of the diagonal of the Oij matrix.
    """
    sum_metric_value_same_clusters = 0
    for cluster_pair in cluster_pairs:
        if(cluster_pair[0] == cluster_pair[1]):
            sum_metric_value_same_clusters += cluster_pairs[
                                        cluster_pair]["metric_value"]
    return sum_metric_value_same_clusters

def compute_expected_metric_value_all_pairs(cluster_pairs):
    """
    Return the expected metric value for each cluster pair.

    This is computed in such a way to preserve
    the total metric value of each cluster.
    The process is: for each cluster pair,
    multiplying together their sums of metric values,
    and then dividing by the total metric values sums of all clusters.
    It is stored in a dictionary field, indexed by the cluster pair.
    """
    expected_metric_values = collections.Counter()
    for cluster_pair in cluster_pairs:
        # This computation corresponds to multiplying together
        # the row (or column) sums corresponding to each cluster
        # in the Oij matrix, and then dividing by the total matrix sum.
        # The denominator in this division is implemented in this way
        # to make use of the existing data structure (cluster_pairs)
        # and is equivalent to multiplying by two the lower triangular
        # + diagonal sum, and then subtracting the sum of the diagonal
        # (remember that Oij is symmetrical)
        expected_metric_values[cluster_pair] = float(
                float(compute_sum_metric_value_cluster(cluster_pair[0]
                                                        , cluster_pairs))
                * float(compute_sum_metric_value_cluster(cluster_pair[1]
                                                        , cluster_pairs))
                / (2.0
                    *float(compute_sum_metric_value_all_clusters(
                                                            cluster_pairs))
                    - float(compute_sum_metric_value_same_clusters(
                                                        cluster_pairs))
                    )
            )
    return expected_metric_values

def random_pctg_sample_of_clones(clones):
    """Return a 30% sample of clones."""
    sample_pctg = 0.3
    sample_size = int(sample_pctg*len(clones))
    clones_sample = random.sample(list(clones.keys()), k=sample_size)
    return clones_sample

def dict_to_array(cluster_pairs, num_clusters, key_arg):
    """
    Return a 2D list with a value for each cluster pair, and their IDs.

    The value for each cluster pair is obtained from the key
    named key_arg.
    The 2D list to be returned is conceptually considered to be
    a 2D array, with the rows and columns corresponding to
    the first and second clusters of the cluster pairs,
    and viceversa (to fill the whole matrix).
    The cell value will be the corresponding item stored in the field
    key_name.
    The list of row/column IDs is in insertion order, which is the same
    that the order in which they appear in a cluster pair
    when the latter are iterated over.
    Obs:
        The latter means that the order of the cluster IDs (idents)
        is not maintained, but it depends on
        when it appears within a cluster pair,
        and when that cluster pair appears in the iteration of
        cluster pairs.
    """
    # Ordered Dictionary that will store the indexes in the array
    # corresponding to each of the clusters
    array_indexes = collections.OrderedDict()
    max_array_index = -1
    # Initialize the array to a (num_clusters x num_clusters) shape
    array_2D = np.array([
                                [0.0 for i in range(num_clusters)]
                            for j in range(num_clusters)])
    for cluster_pair in cluster_pairs:
        for cluster in cluster_pair:
            if(cluster not in array_indexes):
                max_array_index += 1
                print("array_indexes[{}] = {}"
                                            .format(cluster, max_array_index))
                # As array_indexes is an Ordered Dict, the order in
                # which the indexes were inserted (and the dict keys
                # created) will correspond to the indexes themselves
                array_indexes[cluster] = max_array_index
        # Fill the matrix with the corresponding value
        array_2D[array_indexes[cluster_pair[0]]
                ][array_indexes[cluster_pair[1]]] = (
                                        cluster_pairs[cluster_pair][key_arg]
                        )
        # In both orders (to fill the whole matrix)
        array_2D[array_indexes[cluster_pair[1]]
                ][array_indexes[cluster_pair[0]]] = (
                                        cluster_pairs[cluster_pair][key_arg]
                        )
    # The keys in insertion order of array_indexes are the
    # clusterIDs in the order that they will appear in the new 2D array
    keys_in_insertion_order = [key for key in array_indexes]
    return array_2D, keys_in_insertion_order

def print_matrix_onto_csv(matrix, clusters_names_array_order
    , csv_matrix_output_file):
    """
    Write a matrix of values of each cluster pair into a csv file.
    """
    # print_2D_float_list(matrix)
    with open(csv_matrix_output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write the column indexes in the first row
        writer.writerow(["ident_name"] + clusters_names_array_order)
        # Write each row, beginning by its corresponding index
        for cluster_name, row in zip(
                                clusters_names_array_order, matrix):
            writer.writerow([cluster_name] + list(row))
        csvfile.close()

def write_matrix(title, fieldname, cluster_pairs, clusters, clusters_names
    , csv_matrix_output_file):
    """
    Write a matrix of values of each cluster pair into a csv file.

    Also, return the matrix in a 2D list form,
    and a list with the cluster names.
    Parameters:
        title: string that contains the title of the matrix to be
        outputed. It's just for outputing a message in the console.
        fieldname: in the dictionary cluster_pairs, it's the
        name of the field whose values are desired to be written.
        cluster_pairs: dictionary indexed by the cluster pairs,
        among whose attributes should be fieldname.
        clusters: a collections.Counter() object with the
        frequency distribution of the clusters.
        clusters_names: a list that stores, for each clusterID,
        its corresponding name.
        csv_matrix_output_file: string with the file path where
        the matrix is desired to be written.
    """
    # First, transform the data
    # to a structure accepted by the used libraries.

    # Create a 2-D array form of the values in cluster_pairs,
    # and a list that stores the clusters in the order they appear
    # in that array
    matrix, clusters_list_array_order = (
                dict_to_array(cluster_pairs, len(clusters)
                                            , fieldname)
            )
    print(title + ":")
    print_2D_float_list(matrix)
    # Create a list that stores the cluster names, in the order
    # in which they occur in the 2D array (along the columns or rows)
    clusters_names_array_order = [
            clusters_names[cluster] for cluster in clusters_list_array_order]
    print_matrix_onto_csv(matrix, clusters_names_array_order
                            , csv_matrix_output_file)
    return matrix, clusters_names_array_order

def compute_distributions(csv_data_clone_cells, num_iterations
    , cluster_pairs, clusters, clusters_names, clones):
    """
    Return distribution median for each cluster pair.

    Do the num_iterations iterations of random 30% samples of clones
    and subsequent simulation, and return the median metric value of
    each pair of clusters.
    """
    # Dictionary that will store in lists the values
    # obtained over the iterations for each cluster pair
    cluster_pairs_iterations = generate_cluster_pairs(clusters.keys())
    for cluster_pair in cluster_pairs_iterations:
        cluster_pairs_iterations[
                    cluster_pair]["ratio_observed_vs_expected"] = []
    # Dictionary, local to each iteration,
    # that will store the observed metric value of each cluster pair
    cluster_pairs_local_iteration = generate_cluster_pairs(clusters.keys())
    for iteration in range(0, num_iterations):
        # List that will store the data for the current iteration
        csv_data_clone_cells_local_iteration = []
        # List that stores the clones obtained from the random sample
        clones_sample = random_pctg_sample_of_clones(clones)
        # Extract rows whose cloneID belongs to the clones sample
        for csv_data_clone_cell in csv_data_clone_cells:
            if(csv_data_clone_cell["cloneID"] in clones_sample):
                csv_data_clone_cells_local_iteration.append(
                                                        csv_data_clone_cell)        
        # Compute the observed metric value for each cluster pair,
        # according to the data of the current iteration.
        num_cells = compute_num_cells_per_clones_per_cluster(
                                        csv_data_clone_cells_local_iteration)
        cluster_pairs_local_iteration = (
                        compute_observed_metric_value_all_pairs(num_cells
                                    , cluster_pairs_local_iteration, clones)
                )
        # Compute the expected metric value
        # for each cluster pair,
        # according to the data of the current iteration.
        expected_metric_values = (
                        compute_expected_metric_value_all_pairs(
                                                cluster_pairs_local_iteration)
                        )
        # Append the current iteration's results to the lists
        for cluster_pair in cluster_pairs_local_iteration:
            if(expected_metric_values[cluster_pair] > 0.0):
                cluster_pairs_iterations[
                        cluster_pair]["ratio_observed_vs_expected"].append(
                                float(cluster_pairs_local_iteration[
                                        cluster_pair]["metric_value"])
                                / expected_metric_values[cluster_pair]
                )
    # Compute the median of the ratios Oij/Eij
    # observed along the iterations
    for cluster_pair in cluster_pairs_iterations:
        cluster_pairs[cluster_pair]["median_ratio_observed_vs_expected"] = (
                np.median(cluster_pairs_iterations[
                                cluster_pair]["ratio_observed_vs_expected"])
                )
    return cluster_pairs

def print_2D_float_list(l):
    """Print a 2D list of floats in a row-wise format."""
    for row in l:
        for col in row:
            print("{:8.1f}".format(col), end=" ")
        print("")

def plot_clustermap(array_2D, labels_array_order, vmin, vmax):
    """
    Plot clustermap of the values in array_2D.

    An optional step of hiding the upper triangular matrix
    can be taken.
    """
    # print_2D_float_list(array_2D)
    mask = np.zeros_like(array_2D)
    # Set the next assigment to False for plotting
    # the entire matrix, and to True for
    # only the lower triangular matrix with the diagonal
    mask[np.triu_indices_from(mask)] = False
    mask[np.diag_indices_from(mask)] = False
    with sns.axes_style("white"):
        ax = sns.clustermap(array_2D
                    , xticklabels = labels_array_order
                    , yticklabels = labels_array_order
                    , mask=mask, vmin=vmin, vmax=vmax
                    # , annot=True
                    # , square=True
                    ,  cmap="YlOrRd")
        plt.setp(ax.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.show()

def output_lineage_coupling_computations(cluster_pairs, clusters
    , clusters_names, csv_lineage_coupling_scores_matrix_output_file
    , vmin, vmax):
    """
    Output computations related to clusters' lineage couplings.

    The computations of the lineage coupling scores as:
        - a matrix onto a csv
        - and a clustermap in a plot
    """
    # First, transform data to a structure accepted by seaborn methods.
    lineage_coupling_scores_array, clusters_names_array_order = write_matrix(
            "Lineage Coupling Matrix", "median_ratio_observed_vs_expected"
            , cluster_pairs, clusters, clusters_names
            , csv_lineage_coupling_scores_matrix_output_file)
    plot_clustermap(lineage_coupling_scores_array, clusters_names_array_order
                    , vmin, vmax)

if __name__ == "__main__":
    # Obtain the parameters of the program
    (num_iterations, csv_data_url, csv_data_file
        , csv_lineage_coupling_scores_matrix_output_file
        , csv_metric_values_matrix_output_file
        , vmin, vmax) = parse_args()
    # The url is given priority as the data source over the local file
    csv_data_source_str = (csv_data_url
                                if (csv_data_url is not None)
                                else csv_data_file)
    print("Run 'python lineage_coupling_analysis_weinreb_way_5.py -h'"
                                                    " for displaying usage.")
    # Iterable that stores the parsed data
    csv_data_clone_cells = parse_csv_data(csv_data_url, csv_data_file)
    # Counter of cells per clone
    clones = compute_clones_distribution(csv_data_clone_cells)
    # Counter of cells per cluster, and names of each cluster
    clusters, clusters_names = compute_clusters_distribution(
                                    csv_data_clone_cells)
    # Dictionary-valued dictionary indexed by the cluster pairs.
    cluster_pairs = generate_cluster_pairs(clusters.keys())
    # Multi-level counter with different cell counts per level
    num_cells = compute_num_cells_per_clones_per_cluster(csv_data_clone_cells)
    cluster_pairs = compute_observed_metric_value_all_pairs(
                                            num_cells, cluster_pairs, clones)
    print("Writing metric values per pair into csv file '{}'"
            .format(csv_metric_values_matrix_output_file))
    write_matrix("Matrix of metric values per pair"
                , "metric_value", cluster_pairs
                , clusters, clusters_names
                , csv_metric_values_matrix_output_file)
    print("Computing lineage coupling scores for data in '{}'"
        ", with {} iterations, clustermap_minimum {}, and clustermap_maximum {}"
        .format(csv_data_source_str, num_iterations, vmin, vmax))
    cluster_pairs = compute_distributions(csv_data_clone_cells, num_iterations
                        , cluster_pairs, clusters, clusters_names, clones)
    pp.pprint(cluster_pairs)
    output_lineage_coupling_computations(cluster_pairs, clusters
            , clusters_names, csv_lineage_coupling_scores_matrix_output_file
            , vmin, vmax)