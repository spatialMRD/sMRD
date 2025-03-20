import numpy as np
import pandas as pd
from itertools import combinations

DATA_DIR = "data/"
SAMPLE_INFO_PATH = f"{DATA_DIR}/sample_info/ex_sample_info.csv"
CNV_PATH = f"{DATA_DIR}/features/ex_cnv.csv"
SNV_PATH = f"{DATA_DIR}/features/ex_snv.csv"
OUTPUT_DIR = f"{DATA_DIR}/output"


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
    mrd_positive_sample_ids = mrd_positive_sample_ids["sample_id"].tolist()

    # get the pretreatment tumors
    patient_ids = sample_info["patient_id"].unique()
    pretreatment_tumors = {}
    for patient_id in patient_ids:
        tmp = sample_info.loc[sample_info["patient_id"] == patient_id, :]
        tmp = tmp.loc[tmp.sample_type == "Pre-treatment Tumor", :]
        pretx_tumors = tmp.sample_id.tolist()
        pretreatment_tumors[patient_id] = pretx_tumors

    return sample_info, mrd_positive_sample_ids, pretreatment_tumors


# read the cnv file
def load_cnv_data():
    cnv = pd.read_csv(CNV_PATH)
    return cnv


def load_snv_data():
    snv = pd.read_csv(SNV_PATH)
    return snv


def load_features():
    cnv = load_cnv_data()
    snv = load_snv_data()
    # concatenate both dataframes, columnwise
    cnv = cnv.set_index("sample")
    snv = snv.set_index("sample")
    df = cnv.join(snv, how="outer")
    return df


def jaccard_similarity(df):
    """
    Compute the Jaccard similarity between all pairs of rows in a dataframe.

    Args:
        df (pd.DataFrame): Input dataframe where each row represents a sample,
                           and values can be -1, 0, or 1.

    Returns:
        pd.DataFrame: A symmetric dataframe containing
                      Jaccard similarities between rows.
    """
    # Convert dataframe to numpy array for faster computations
    data = df.values

    # Get the number of rows
    n_rows = data.shape[0]

    # Initialize an empty similarity matrix
    similarity_matrix = np.zeros((n_rows, n_rows))

    # Compute pairwise Jaccard similarities
    for i, j in combinations(range(n_rows), 2):
        intersection = np.sum((data[i] != 0) & (data[j] != 0) & (data[i] == data[j]))
        union = np.sum((data[i] != 0) | (data[j] != 0))
        similarity = intersection / union if union != 0 else 0

        # Update similarity matrix
        similarity_matrix[i, j] = similarity
        similarity_matrix[j, i] = similarity

    # Set diagonal elements to 1 (self-similarity)
    np.fill_diagonal(similarity_matrix, 1.0)

    # Convert to a DataFrame for easier interpretation
    similarity_df = pd.DataFrame(similarity_matrix, index=df.index, columns=df.index)

    return similarity_df


def get_jaccard_stats(df, patient_id, mrd_positive_sample_ids, pretreatment_tumors):
    # print stats on pretreatment tumors
    samples = [x for x in df.index if patient_id in x]
    samples = [x for x in samples if x in mrd_positive_sample_ids]
    patient_df = df.loc[samples]
    patient_df = patient_df.fillna(0)
    similarity = jaccard_similarity(patient_df)
    # drop diagonals from similarity
    similarity = similarity.replace(1, np.nan)
    # flatten
    tmp = similarity.values.flatten()
    # drop nans
    tmp = tmp[~np.isnan(tmp)]
    stats = {
        "patient_id": patient_id,
        "mean": tmp.mean(),
        "std": tmp.std(),
        "min": tmp.min(),
        "max": tmp.max(),
        "median": np.median(tmp),
    }
    return stats


def main():
    sample_info, mrd_positive_sample_ids, pretreatment_tumors = load_sample_info()
    features = load_features()

    patient_ids = sample_info["patient_id"].unique()

    patient_stats = []
    for patient_id in patient_ids:
        patient_stats.append(
            get_jaccard_stats(
                features, patient_id, mrd_positive_sample_ids, pretreatment_tumors
            )
        )

    jaccard_df = pd.DataFrame(patient_stats)
    # add metadata to jaccard_df
    patient_dict = dict(zip(sample_info["patient_id"], sample_info["progression"]))
    jaccard_df["progression"] = jaccard_df["patient_id"].map(patient_dict)

    jaccard_df.to_csv(f"{OUTPUT_DIR}/jaccard_stats.csv", index=False)


if __name__ == "__main__":
    main()
