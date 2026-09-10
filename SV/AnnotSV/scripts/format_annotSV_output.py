import pandas as pd
import numpy as np

def filter_samples_annotSV_output(annotSV_output_path, samples_file_path, output_path):
    """
    Filter AnnotSV output to include samples in the samples file and keep only the columns specified
    """

    annotSV_all = pd.read_csv(annotSV_output_path, sep="\t")

    annotSV_cols = ['SV_chrom', 'SV_start', 'SV_end', 'Samples_ID', 'id', 'SV_type', 'Annotation_mode', 'CytoBand', 'Gene_name', 'AnnotSV_ID']
    annotSV_all = annotSV_all[annotSV_cols]

    samples = pd.read_csv(samples_file_path)["tumor_normal_pair"].tolist()

    # comparison stats
    print(
        "Number of samples in annotSV_output not in samples to include:",
        annotSV_all.loc[~annotSV_all["Samples_ID"].isin(samples), "Samples_ID"].nunique()
    )
    print(
        "Number of samples to include not in annotSV_output:",
        len(set(samples) - set(annotSV_all["Samples_ID"]))
    )
    # print the SVs that are in samples but not in annotSV_output
    missing_samples = sorted(set(samples) - set(annotSV_all["Samples_ID"].dropna().astype(str)))
    print("\n".join(missing_samples))

    # filter and save
    annotSV_output = annotSV_all[annotSV_all["Samples_ID"].isin(samples)]
    annotSV_output.to_csv(output_path, sep="\t", index=False)

if __name__ == "__main__":
    annotSV_output_path = "../results/sv_calls_v2.annotated.tsv"
    analysis_samples_path = "../../SV_sample_data/data/final_analysis_samples.csv"
    final_samples_path = "../../SV_sample_data/data/final_samples.csv"
    filter_samples_annotSV_output(annotSV_output_path, final_samples_path, "../results/final_samples_annotSV_output.tsv")
    filter_samples_annotSV_output(annotSV_output_path, analysis_samples_path, "../results/analysis_samples_annotSV_output.tsv")