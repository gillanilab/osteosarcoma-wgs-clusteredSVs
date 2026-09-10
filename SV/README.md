SV_sample_data
==============

### Workflow
1. `SV_sample_data/sv.ipynb` - analyze sv calling, save summary files
2. `SV_sample_data/sv_landscape_plot.ipynb` - (Figure 2C; Extended Data Figure 1B)
3. `SV_sample_data/plot_svs.R` - plot rearrangement profiles using ReConPlot (Figure 4C, 4D; Extended Data Figures 1C, 3A)
4. `SV_sample_data/data_prep_for_oncoprint.ipynb`
5. `SV_sample_data/oncoprint.R` - (Figure 1B)


AnnotSV (SV Annotation)
=======================

### Requirements
- AnnotSV v3.4.6 installed locally at: SV/tools/AnnotSV (https://github.com/lgmgeo/AnnotSV)

### Workflow
1. `AnnotSV/scripts/h5_to_bed.py` — convert SV calls from HDF5 to BED
2. `AnnotSV/scripts/run_annotsv.sh` — run AnnotSV
3. `AnnotSV/scripts/format_annotSV_output.py` - filter and format AnnotSV output

### Notes
- For SV type, breakpoints from long_range and inter_chr SVs are entered separately and labeled as `BND`
- AnnotSV typically takes ~15 minutes to run
 

SRB_Analysis (Significantly Recurrent Bins Analysis)
====================================================

### Workflow
1. `SRB_analysis/scripts/fishHook_installation.R`
2. `SRB_analysis/notebooks/bins.ipynb` - divide breakpoints across genomic bins
3. `SRB_analysis/notebooks/covariates.ipynb` - compute per-bin genomic covariates
4. `SRB_analysis/notebooks/background_model.ipynb` - merge breakpoint and covariate data, prepare inputs for fishhook
5. `SRB_analysis/scripts/run_fishhook.R`
6. `SRB_analysis/notebooks/interpret_fishhook.ipynb` - interpret & visualize results (Figure 2A)


Shatterseek (Chromothripsis)
============================

### Workflow
1. `Shatterseek/scripts/shatterseek_installation.R`
2. `Shatterseek/notebooks/data_prep.ipynb` — prepare the SV and CNV files for Shatterseek
3. `Shatterseek/scripts/run_shatterseek.R` - run Shatterseek
4. `Shatterseek/notebooks/intepret_shatterseek.ipynb` - explore Shatterseek output, classify CGR and chromothripsis events, annotate events with genes
5. `Shatterseek/scripts/plot_shatterseek.R` — plots specified events


SV_signature_analysis
=====================

### Requirements
- Signature tools installed locally at: SV/tools/signature.tools.lib (https://github.com/Nik-Zainal-Group/signature.tools.lib)

    Clustered_SVs
    =============

    ### Workflow
    1. `SV_signature_analysis/Clustered_SVs/notebooks/data_prep.ipynb` - prepare bedpe as input for signature tools package
    2. `SV_signature_analysis/Clustered_SVs/scripts/clustered.R` - classify SVs as clustered/not clustered & find clustering regions
    3. `SV_signature_analysis/Clustered_SVs/notebooks/identify_hotspots.ipynb` - identify the clustered SV hotspot regions
    4. `SV_signature_analysis/Clustered_SVs/notebooks/linking_clustered_svs.ipynb` - find linked clustered SVs & merge linked intra-chrom clustered SV regions

    ### Other Analyses
    - `SV_signature_analysis/Clustered_SVs/scripts/SV_counts_comparison` - comparing SV counts in primary vs metastatic samples (Figure 2D)
    - `SV_signature_analysis/Clustered_SVs/notebooks/classify_clustered_sv.ipynb` - find overlap between clustered SVs, shatterseek CGR, and AA ecDNA (Figure 3B)
    - `SV_signature_analysis/Clustered_SVs/notebooks/lta_analysis.ipynb` - (Figure 3E)
    - `SV_signature_analysis/Clustered_SVs/notebooks/genomic_feature_enrichment.ipynb` - (Extended Data Figure 1E)
    - `SV_signature_analysis/Clustered_SVs/scripts/make_karyoploter.R` - create clustered SV hotspot plot for manuscript (Figure 2E, Extended Data Figure 1D)
    - `SV_signature_analysis/Clustered_SVs/scripts/purity_comparison.R` - (Extended Data Figure 1H)
    - `SV_signature_analysis/Clustered_SVs/scripts/TSG_inactivation_clustered_SVs.R` - associations between TSG disruption and clustered SV events (Figure 3C)
    - `SV_signature_analysis/Clustered_SVs/scripts/clustered_SV_pairwise.R` - pairwise associations between clustered SV events(Figure 3D)

    ### Notes
    - All rearrangements <1kb are ignored by signature.tools.lib package

    TAD_BA_SVs
    ==========
    ### Workflow
    1. `SV_signature_analysis/TAD_BA_SVs/notebooks/tad_disruption.ipynb` - conduct TAD BA-SV analysis
    2. `SV_signature_analysis/TAD_BA_SVs/notebooks/tad_basv_figure.ipynb` - (Extended Data Figure 1F)


CNAs_Amplicon_Architect
=======================

### Analyses
1. run and analyze results of Amplicon Architect (ecDNA)
    - `CNAs_Amplicon_Architect/scripts/amplicon_architect.wdl` - Run workflow on Terra cloud computing platform
    - `CNAs_Amplicon_Architect/notebooks/analysis_aa_results.ipynb`
    - `CNAs_Amplicon_Architect/notebooks/ecDNA_enrichment_analysis.ipynb` - test enrichment of ecDNA at hotspot loci
2. create circos-style plots
    - `CNAs_Amplicon_Architect/notebooks/co-amplification_plots.ipynb` - (Figure 3E, 3F)
    - `CNAs_Amplicon_Architect/notebooks/cnv_circos_plots.ipynb` - (Extended Data Figure 2G)
3. visualize CNAs associated with hotspot clustered SVs
    - `CNAs_Amplicon_Architect/notebooks/hotspot_amplifications.ipynb` - (Figure 3A)


Clonal_analysis
===============

### Workflow
1. `SVClone/scripts/multi_tumor_classification.py` - classification of SVs as truncal, shared, or private
2. `SVClone/notebooks/truncal_sv_analysis.ipynb` - generate input & figure for analysis of % clustered SVs that are truncal across hotspot and non-hotspot regions (Figure 5A)
3. `SVClone/scripts/glmm_truncal_sv_analysis.R` - test whether clustered SVs in hotspot regions are more likely to be truncal using GLMM
4. For SVClone, prepare input files as described on https://github.com/mcmero/SVclone (Copy number inputs were provided in ASCAT format) 
5. `SVClone/scripts/SV_Clone.wdl` (Make sure that path to `svclone_config.ini` is provided)
6. `SVClone/scripts/SV_Clone_postassign.wdl` - get refined clustering based on SNV clustering before. 
7. `SVClone/notebooks/SV_subclonality.ipynb` - find subclonality of hotspot clustered SVs (Figure 5B)
8. `SVClone/notebooks/distinct_clones.ipynb` -  find number of distinct CCF clusters (clones) for hotspot events (Figure 5C)
9. `SVClone/notebooks/KDE_plots.ipynb` - (Figure 5C)


