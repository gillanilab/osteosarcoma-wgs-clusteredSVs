import pandas as pd
import numpy as np
import h5py

def convert_h5_to_bed(h5_file, output_bed):
    df = pd.read_hdf(h5_file).copy()

    df["sv_id"] = np.arange(1, len(df) + 1)

    bed_df = pd.DataFrame()
    bed_df["chrom"] = df["chr1"].astype(str)
    bed_df["start"] = df["pos1"].astype(int)
    bed_df["end"] = bed_df["start"] + 1
    bed_df["name"] = df["name"].astype(str)     # sample
    bed_df["id"] = df["sv_id"].astype(str)      # SV id (pairs breakends)
    bed_df["svtype"] = "BND"

    bed_df.loc[df["class"] == "deletion", "svtype"] = "DEL"
    bed_df.loc[df["class"] == "tandem_dup", "svtype"] = "DUP"
    bed_df.loc[df["class"] == "inversion", "svtype"] = "INV"

    # inter_chr and long_range become BNDs (if we convert long_range to 1 entry of DEL/DUP/INV, the program crashes)
    bed_df.loc[df["class"].isin(["inter_chr", "long_range"]), "svtype"] = "BND"

    # actual end for non-BND
    non_bnd = bed_df["svtype"] != "BND"
    bed_df.loc[non_bnd, "end"] = df.loc[non_bnd, "pos2"].astype(int).values
    s = bed_df.loc[non_bnd, "start"].values
    e = bed_df.loc[non_bnd, "end"].values
    bed_df.loc[non_bnd, "start"] = np.minimum(s, e)
    bed_df.loc[non_bnd, "end"] = np.maximum(s, e)

    # keep end as start + 1 for "BND", add separate entry for partner breakpoint (same name+id)
    bnd = bed_df["svtype"] == "BND"
    if bnd.any():
        partner = bed_df.loc[bnd].copy()
        partner["chrom"] = df.loc[bnd, "chr2"].astype(str).values
        partner["start"] = df.loc[bnd, "pos2"].astype(int).values
        partner["end"] = partner["start"] + 1
        bed_df = pd.concat([bed_df, partner], ignore_index=True)

    bed_df["chrom"] = bed_df["chrom"].replace({"23": "X", "24": "Y"})
    bed_df = bed_df[["chrom","start","end","name","id","svtype"]].sort_values(["chrom","start","end"])

    with open(output_bed, "w") as f:
        f.write("#chrom\tstart\tend\tname\tid\tsvtype\n")
    bed_df.to_csv(output_bed, sep="\t", index=False, header=False, mode="a")

    return bed_df

if __name__ == "__main__":
    h5_file = "../../SV_sample_data/analysis_public_os_data_merged_consensus_svs_v2.h5"
    output_bed = "../data/sv_calls_v2.bed"

    bed_df = convert_h5_to_bed(h5_file, output_bed)