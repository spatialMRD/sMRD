# Compute Spatial and Genomic Correlation
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob

DATA_DIR = "data/"
SAMPLE_INFO_PATH = f"{DATA_DIR}/sample_info/ex_sample_info.csv"
PHYLO_PATH = f"{DATA_DIR}/output/phylo"

def load_sample_info() -> pd.DataFrame:
    sample_info = pd.read_csv(SAMPLE_INFO_PATH)
    # drop prows where progression is NaN
    sample_info = sample_info.dropna(subset=["progression"])
    # set progression col as int
    sample_info["progression"] = sample_info["progression"].astype(int)
    # select only those to include
    sample_info = sample_info.loc[sample_info["include"] == 1]
    # get the MRD positive blocks as a list
    mrd_positive_sample_ids = sample_info.loc[sample_info["detected"] == 1]
    mrd_positive_sample_ids = mrd_positive_sample_ids['sample_id'].tolist()

    # get the pretreatment tumors
    patient_ids = sample_info["patient_id"].unique()
    pretreatment_tumors = {}
    for patient_id in patient_ids:
        tmp = sample_info.loc[sample_info["patient_id"] == patient_id, :]
        tmp = tmp.loc[tmp.sample_type == "Pre-treatment Tumor", :]
        pretx_tumors = tmp.sample_id.tolist()
        pretreatment_tumors[patient_id] = pretx_tumors

    return sample_info, mrd_positive_sample_ids, pretreatment_tumors


def main():
    sample_info, _, pretreatment_tumors = load_sample_info()

    dist_mats = glob(f"{PHYLO_PATH}/*/*_phylo_distances.csv")

    dist = {}
    for fn in dist_mats:
        patient_id = os.path.basename(fn).split("_")[0]
        dist_mat = pd.read_csv(fn, index_col=0)
        dist[patient_id] = dist_mat

    patient_ids = sample_info["patient_id"].unique()

    corrs = []
    for patient_id in patient_ids:
        print(patient_id)
        pretx_tumor = pretreatment_tumors[patient_id][0]
        tmp = dist[patient_id] 
        # get mrd detected for this patient
        mrd_detected = sample_info.loc[sample_info['patient_id'] == patient_id, :]
        mrd_detected = mrd_detected.loc[mrd_detected.detected == 1, :].sample_id.tolist()

        # intersect mrd_detected with tmp.columns
        mrd_detected = list(set(mrd_detected).intersection(tmp.columns))

        print("   Available before MRD filter: ", len(tmp.columns))
        tmp = tmp.loc[mrd_detected, mrd_detected].to_dict()[pretx_tumor]
        print("   Available after MRD filter: ", len(tmp))
        # tmp as dataframe, cols sample_id and cnv_dist
        tmp = pd.DataFrame(tmp.items(), columns=["sample_id", "cnv_dist"])

        # join with spatial
        df = tmp.merge(sample_info, on="sample_id", how="left")

        # drop nan
        print("   Available before dropna: ", len(df))
        df = df.dropna(subset=["x", "y", "z"])
        print("   Available after dropna: ", len(df))
        df = df.sort_values('cnv_dist')

        if len(df) <= 3:
            print("   Skipping patient_id: ", patient_id)
            continue

        # get closest sample
        closest_coords = df.iloc[1, :][['x', 'y', 'z']].values

        # compute the euclidean distance between all samples and the closest sample
        df['euclidean_dist'] = np.sqrt((df['x'] - closest_coords[0])**2 + (df['y'] - closest_coords[1])**2 + (df['z'] - closest_coords[2])**2)

        # get the correlation between euclidean_dist and cnv_dist
        corr = float(df[['euclidean_dist', 'cnv_dist']].corr().iloc[0,1])
        corrs.append({
            "patient_id": patient_id,
            "corr": corr
        })

    corrs = pd.DataFrame(corrs)
    patient_id_to_progression = dict(zip(sample_info["patient_id"], sample_info["progression"]))
    corrs['progression'] = corrs['patient_id'].map(patient_id_to_progression)

    corrs.to_csv(f"{DATA_DIR}/output/spatial_correlation.csv", index=False)

if __name__ == "__main__":
    main()