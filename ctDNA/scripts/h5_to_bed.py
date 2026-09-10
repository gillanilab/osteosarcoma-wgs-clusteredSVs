import pandas as pd
import numpy as np

def convert_h5_to_bed(sv_file, output_bed):
    df = pd.read_csv(sv_file).copy()

    df["sv_id"] = np.arange(1, len(df) + 1)

    bed_df = pd.DataFrame()
    bed_df["chrom"] = df["chr1"].astype(str)
    bed_df["start"] = df["pos1"].astype(int) - 1
    bed_df["end"] = bed_df["start"] + 1
    bed_df["name"] = df["sample"].astype(str)     # sample
    bed_df["id"] = df["sv_id"].astype(str)      # SV id (pairs breakends)

    bed_df.loc[df["svtype"] == "DEL", "svtype"] = "DEL"
    bed_df.loc[df["svtype"] == "DUP/INS", "svtype"] = "DUP"
    bed_df.loc[df["svtype"] == "INV", "svtype"] = "INV"

    # translocations become BNDs
    bed_df.loc[df["svtype"].isin(["TRA"]), "svtype"] = "BND"

    # actual end for non-BND
    non_bnd = bed_df["svtype"] != "BND"
    p1 = df["pos1"].astype(int)
    p2 = df["pos2"].astype(int)
    bed_df.loc[non_bnd, "start"] = np.minimum(p1[non_bnd], p2[non_bnd]) - 1
    bed_df.loc[non_bnd, "end"]   = np.maximum(p1[non_bnd], p2[non_bnd])

    # keep end as start + 1 for "BND", add separate entry for partner breakpoint (same name+id)
    bnd = bed_df["svtype"] == "BND"
    if bnd.any():
        partner = bed_df.loc[bnd].copy()
        partner["chrom"] = df.loc[bnd, "chr2"].astype(str).values
        partner["start"] = df.loc[bnd, "pos2"].astype(int).values - 1
        partner["end"] = partner["start"] + 1
        bed_df = pd.concat([bed_df, partner], ignore_index=True)

    bed_df["chrom"] = bed_df["chrom"].replace({"23": "X", "24": "Y"})
    bed_df = bed_df[["chrom","start","end","name","id","svtype"]].sort_values(["chrom","start","end"])

    with open(output_bed, "w") as f:
        f.write("#chrom\tstart\tend\tname\tid\tsvtype\n")
    bed_df.to_csv(output_bed, sep="\t", index=False, header=False, mode="a")

    return bed_df

if __name__ == "__main__":
    sv_file = "../data/leopard_cohort/final_consensus_sv_hg38_filtered.csv"
    output_bed = "../data/leopard_cohort/annotsv_input.bed"

    bed_df = convert_h5_to_bed(sv_file, output_bed)