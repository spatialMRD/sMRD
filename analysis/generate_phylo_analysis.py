import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skbio import DistanceMatrix
from skbio.tree import nj
from tqdm import tqdm
from Bio import Phylo

DATA_DIR = "data/"
CNV_PATH = f"{DATA_DIR}/features/ex_cnv.csv"
SNV_PATH = f"{DATA_DIR}/features/ex_snv.csv"
OUTPUT_DIR = f"{DATA_DIR}/output/phylo"

def plot_tree(tree, output_file):
    try:
        # Convert the skbio TreeNode to Newick format
        newick_str = tree.write(file="tmp.tree", format="newick")
        
        # Read the tree using Bio.Phylo
        bio_tree = Phylo.read("tmp.tree", "newick")
        
        # Plot the tree
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(bio_tree, do_show=False, axes=ax)
        
        # Save to file
        plt.savefig(output_file, format="png")
        plt.close(fig)
        print(f"Tree saved to {output_file}")
    except Exception as e:
        print(f"Failed to plot tree: {e}")

def single_linkage_clustering(dist_matrix: np.ndarray, root_idx: int) -> np.ndarray:
    """
    Perform single-linkage clustering starting from a specified root node.

    Args:
        dist_matrix (np.ndarray): Square distance matrix
        root_idx (int): Index of the root node to start clustering from

    Returns:
        np.ndarray: Tree distance matrix representing single-linkage clustering results
    """
    # Initialize sets
    connected_nodes = [root_idx]
    remaining_nodes = list(range(dist_matrix.shape[0]))
    remaining_nodes.remove(root_idx)

    # Initialize a result distance matrix to store tree connections, fill with np.inf
    tree_dist_matrix = np.full((dist_matrix.shape[0], dist_matrix.shape[0]), np.inf)

    # Perform the iterative process
    while remaining_nodes:
        closest_node = None
        closest_dist = np.inf
        parent_node = None

        # Find the closest node to the tree based on asymmetrical distances
        for pnode in connected_nodes:
            for cnode in remaining_nodes:
                dist = dist_matrix[pnode, cnode]
                if dist < closest_dist:
                    closest_dist = dist
                    closest_node = cnode
                    parent_node = pnode

        # Add the closest node to the tree
        connected_nodes.append(closest_node)
        remaining_nodes.remove(closest_node)

        # Update the tree distance matrix (asymmetrical)
        tree_dist_matrix[parent_node, closest_node] = closest_dist
        tree_dist_matrix[closest_node, parent_node] = closest_dist

    return tree_dist_matrix

def hamming_distance(sample1, sample2):
    return sum(abs(a - b) for a, b in zip(sample1, sample2))

def root_tree(tree, root_sample):
    root_node = tree.find(root_sample)
    if root_node is None:
        raise ValueError(f"Sample {root_sample} not found in the tree.")
    tree = tree.root_at(root_node, reset=True, branch_attrs=[], root_name=root_sample)
    return tree

def calculate_phylogenetic_distances(tree):
    tip_names = [tip.name for tip in tree.tips()]
    distances = pd.DataFrame(index=tip_names, columns=tip_names, dtype=float)
    for name1 in tip_names:
        for name2 in tip_names:
            if name1 == name2:
                distances.loc[name1, name2] = 0.0
            else:
                distances.loc[name1, name2] = tree.find(name1).distance(tree.find(name2))
    return distances

def create_trees(df, patient_id):
    samples = df.index
    num_samples = len(samples)

    # create pairwise distance matrix
    distance_matrix = []
    for i in range(num_samples):
        row = []
        for j in range(num_samples):
            if i == j: 
                row.append(0.0)
            else:
                row.append(hamming_distance(df.iloc[i], df.iloc[j]))
        distance_matrix.append(row)

    dm = DistanceMatrix(distance_matrix, ids=samples)
    tree = nj(dm)

    tumor_samples = [sample for sample in samples if "-T" in sample]

    phylogenetic_distances = calculate_phylogenetic_distances(tree)

    if len(tumor_samples) != 0:
        root_idx = list(phylogenetic_distances.index).index(tumor_samples[0])
        cluster = single_linkage_clustering(phylogenetic_distances.values, root_idx)
        cluster = pd.DataFrame(cluster, index=phylogenetic_distances.index, columns=phylogenetic_distances.columns)

        # root the tree
        rooted_tree = root_tree(tree, tumor_samples[0])
        patient_dir = os.path.join(OUTPUT_DIR, patient_id)
        tree_plot_file = os.path.join(patient_dir, "tree.png")
        plot_tree(rooted_tree, tree_plot_file)
    else:
        print("MISSING ROOT TUMOR!")
        root_idx = 0
        cluster = single_linkage_clustering(phylogenetic_distances.values, root_idx)
        cluster = pd.DataFrame(cluster, index=phylogenetic_distances.index, columns=phylogenetic_distances.columns)

    return tree, phylogenetic_distances, cluster

def load_cnv_data():
    cnv = pd.read_csv(CNV_PATH)
    # prepend column for patient_id
    cnv["patient_id"] = cnv["sample"].str.split("-", expand=True)[0]
    # set patient_id as first column
    cnv = cnv[["patient_id"] + [col for col in cnv.columns if col != "patient_id"]]
    return cnv

def load_snv_data():
    snv = pd.read_csv(SNV_PATH)
    # preend column for patient_id
    snv['patient_id'] = snv['sample'].str.split('-', expand=True)[0]
    # set patient_id as first column
    snv = snv[["patient_id"] + [col for col in snv.columns if col != "patient_id"]]
    return snv

def select_patient_data(cnv, patient_id):
    df = cnv.loc[cnv["patient_id"] == patient_id, :]
    # drop patient_id column, make sample column the index
    df = df.drop(columns="patient_id")
    df = df.set_index("sample")
    # compute the number of samples where each variant != 0
    num_samples = df.astype(bool).sum(axis=0)
    # get where num_samples > 0
    num_samples = num_samples[num_samples > 0]
    valid_vars = num_samples.index
    df = df.loc[:, valid_vars]
    return df

def main():
    cnv = load_cnv_data()
    snv = load_snv_data()
    unique_patients = cnv["patient_id"].unique()

    # concatenate both dataframes, columnwise
    cnv = cnv.set_index("sample")
    snv = snv.set_index("sample")
    # drop patient_id from snv
    snv = snv.drop(columns="patient_id")
    cnv = cnv.join(snv, how="outer")
    # reset sample index
    cnv = cnv.reset_index()

    for patient_id in tqdm(unique_patients):
        print(f"Processing patient {patient_id}")
        patient_df = select_patient_data(cnv, patient_id)
        if len(patient_df) < 3:
            print(f"Skipping patient {patient_id} with only {len(patient_df)} samples")
            continue
        tree, phylogenetic_distances, cluster = create_trees(patient_df, patient_id)

        # save output to file
        patient_output = os.path.join(OUTPUT_DIR, patient_id)
        if not os.path.exists(patient_output):
            os.makedirs(patient_output)
        tree.write(f"{patient_output}/{patient_id}_tree.nwk")
        phylogenetic_distances.to_csv(f"{patient_output}/{patient_id}_phylo_distances.csv")
        cluster.to_csv(f"{patient_output}/{patient_id}_cluster.csv")

if __name__ == "__main__":
    main()