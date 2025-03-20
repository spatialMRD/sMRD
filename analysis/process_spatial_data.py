import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm

DATA_DIR = "data/"
SAMPLE_INFO_PATH = os.path.join(DATA_DIR, "sample_info/ex_sample_info.csv")
DIST_MAT_DIR = os.path.join(DATA_DIR, "output/phylo")

def load_sample_info():
    sample_info = pd.read_csv(SAMPLE_INFO_PATH)
    # drop prows where progression is NaN
    sample_info = sample_info.dropna(subset=["progression"])
    # set progression col as int
    sample_info["progression"] = sample_info["progression"].astype(int)
    # select only those to include
    sample_info = sample_info.loc[sample_info["include"] == 1]
    # get a dictionary of patient_id to progression status
    patient_dict = dict(zip(sample_info["patient_id"], sample_info["progression"]))
    # get the MRD positive blocks as a list
    mrd_positive_sample_ids = sample_info.loc[sample_info["detected"] == 1]
    mrd_positive_sample_ids = mrd_positive_sample_ids['sample_id'].tolist()
    # get a dictionary of patient_id to progression status
    patient_id_to_progression = dict(zip(sample_info['patient_id'], sample_info['progression']))
    # get the pretreatment tumors
    # get pretreatment tumors
    patient_ids = sample_info['patient_id'].unique()
    pretreatment_tumors = {}
    for patient_id in patient_ids:
        tmp = sample_info.loc[sample_info['patient_id'] == patient_id, :]
        tmp = tmp.loc[tmp.sample_type == "Pre-treatment Tumor", :]
        pretx_tumors = tmp.sample_id.tolist()
        # if len > 1, get the one with "L"
        if len(pretx_tumors) > 1:
            for tumor in pretx_tumors:
                if "L" in tumor:
                    pretx_tumors = [tumor]
                    break
        pretreatment_tumors[patient_id] = pretx_tumors

    return sample_info, patient_dict, mrd_positive_sample_ids, patient_id_to_progression, pretreatment_tumors

def load_dist_mats():
    # Now load the phylo output
    dist_mats = glob(os.path.join(DIST_MAT_DIR, "*/*_phylo_distances.csv"))
    # load each of these dist mats
    dist = {}
    for fn in dist_mats:
        # get the patient_id from the dist mat
        patient_id = os.path.basename(fn).split("_")[0]
        # load the dist mat
        dist_mat = pd.read_csv(fn, index_col=0)
        # add to the dict
        dist[patient_id] = dist_mat
    return dist

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


def get_spatial_plot(patient_id, patient_dist, pretx_tumor, sample_info):
    backup_patient_dist = patient_dist.copy()
    mrd_detected = sample_info.loc[sample_info["patient_id"] == patient_id, :]
    mrd_detected = mrd_detected.loc[mrd_detected.detected == 1, :].sample_id.tolist()

    # intersect mrd_detected and patient_dist.columns
    mrd_detected = list(set(mrd_detected).intersection(patient_dist.columns))
    patient_dist = patient_dist.loc[mrd_detected, mrd_detected].to_dict()[pretx_tumor]
    patient_dist = pd.DataFrame(patient_dist.items(), columns=["sample_id", "genomic_dist"])

    # join with spatial
    df = patient_dist.merge(sample_info, on="sample_id", how="left")

    # drop nan
    df = df.dropna(subset=["x", "y", "z"])
    df = df.sort_values('genomic_dist')

    # get closest sample
    closest_sample = df.iloc[1, :].sample_id

    # get list of samples we still have
    samples_with_mrd_and_coords = df.sample_id.tolist()

    # filter distance matrix to this
    all_dist = backup_patient_dist.loc[samples_with_mrd_and_coords, samples_with_mrd_and_coords]

    # perform single linkage clustering
    closest_sample_idx = samples_with_mrd_and_coords.index(closest_sample)

    clustering = single_linkage_clustering(all_dist.values, closest_sample_idx)
    clustering = pd.DataFrame(clustering, index=samples_with_mrd_and_coords, columns=samples_with_mrd_and_coords)

    # get the connections
    connections = []
    for i in range(clustering.shape[0]):
        for j in range(i + 1, clustering.shape[0]):
            if clustering.iloc[i, j] != np.inf:
                connections.append((clustering.index[i], clustering.columns[j], clustering.iloc[i, j]))

    df_plot = df.loc[:, ["sample_id", "x", "y", "z", "genomic_dist"]]

    return df_plot, connections

def add_pcs_and_scale(df):
    scaler = MinMaxScaler()
    df["genomic_dist_scaled"] = scaler.fit_transform(df[["genomic_dist"]])

    # Extract x, y, z coordinates
    coordinates = df[["x", "y", "z"]]

    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(coordinates)

    # Standard scale PC1 and PC2
    standard_scaler = StandardScaler()
    pca_result_scaled = standard_scaler.fit_transform(pca_result)
    df["PC1"] = pca_result_scaled[:, 0]
    df["PC2"] = pca_result_scaled[:, 1]

    # Shift the data so the min genomic_dist sample is centered on (0, 0)
    min_dist_index = df["genomic_dist_scaled"].idxmin()
    df["PC1"] -= df.loc[min_dist_index, "PC1"]
    df["PC2"] -= df.loc[min_dist_index, "PC2"]

    return df

def main():
    # load the sample info
    sample_info, patient_dict, mrd_positive_sample_ids, patient_id_to_progression, pretreatment_tumors = load_sample_info()
    # load the dist mats
    dist = load_dist_mats()
    # get the patient ids
    patient_ids = list(dist.keys())

    processed_data = {}
    for patient_id in patient_ids:
        try:
            pretx_tumor = pretreatment_tumors[patient_id][0]
            tmp_dist = dist[patient_id]

            df_plot, connections = get_spatial_plot(patient_id, tmp_dist, pretx_tumor, sample_info)

            df_plot = add_pcs_and_scale(df_plot)

            processed_data[patient_id] = (df_plot, connections)
        except:
            print(f"Failed for {patient_id}")

    # get progressors and non progressors
    progs = sample_info.loc[sample_info.progression == 1, :].patient_id.unique()
    nonprogs = sample_info.loc[sample_info.progression == 0, :].patient_id.unique()

    # filter to only those in processed_data
    progs = list(set(progs).intersection(processed_data.keys()))
    nonprogs = list(set(nonprogs).intersection(processed_data.keys()))

    prog_dfs = []
    prog_connections = []
    for prog in progs:
        df_plot, connections = processed_data[prog]
        prog_dfs.append(df_plot)
        prog_connections.append(connections)

    prog_dfs = pd.concat(prog_dfs)
    prog_connections = [item for sublist in prog_connections for item in sublist]

    nonprog_dfs = []
    nonprog_connections = []
    for nonprog in nonprogs:
        df_plot, connections = processed_data[nonprog]
        nonprog_dfs.append(df_plot)
        nonprog_connections.append(connections)

    nonprog_dfs = pd.concat(nonprog_dfs)
    nonprog_connections = [item for sublist in nonprog_connections for item in sublist]

    prog_dfs["progression"] = 1
    nonprog_dfs["progression"] = 0
    all_dfs = pd.concat([prog_dfs, nonprog_dfs])
    all_connections = prog_connections + nonprog_connections

    # save to file
    all_dfs.to_csv(os.path.join(DATA_DIR, "output/all_spatial_dfs.csv"), index=False)
    with open(os.path.join(DATA_DIR, "output/all_spatial_connections.txt"), "w") as f:
        for conn in all_connections:
            f.write(f"{conn[0]},{conn[1]},{conn[2]}\n")

if __name__ == "__main__":
    main()
