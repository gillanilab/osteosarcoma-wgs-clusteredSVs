ctDNA
=====

### Workflow
1. Post-processing of SV calls (from SvABA), generate input files for clustered SV calling
    - `ctDNA/notebooks/cmi_cohort/cmi_sv_summary.ipynb` - ctDNA SVs
    - `ctDNA/notebooks/leopard_cohort/leopard_sv_summary.ipynb` - ctDNA and matched tumor SVs
2. Clustered SV calling
    - `SV/SV_signature_analysis/Clustered_SVs/scripts/clustered.R`
    - `ctDNA/notebooks/leopard_cohort/merge_clustered_svs.ipynb` - merged connected clustered SV regions into clustered SV events
3. Annotate SVs with genes using AnnotSV to identify TSG-disrupting SVs
    -  Install AnnotSV v3.4.6 locally at: SV/tools/AnnotSV (https://github.com/lgmgeo/AnnotSV)
    - `ctDNA/scripts/h5_to_bed.py` - generate input file for AnnotSV
    - `ctDNA/scripts/run_annotSV.sh` 
4. Visualize SVs
    - `ctDNA/notebooks/circos_plot.ipynb` (Figure 5E, 5F, 5G, 5H, 5J; Extended Data Figure 3B, 3C, 3D)

