SNV Analysis
============

### Workflow
1. `SNV/Signature_Analyzer/notebooks/de_novo_signature_analysis.ipynb`
2. `SNV/Signature_Analyzer/notebooks/interpret_signatures.ipynb` - identify artifact signatures and filter SNVs with high artifact contribution
3.  Annotate SNVs through Genome Nexus Annotation Pipeline (https://github.com/genome-nexus/genome-nexus-annotation-pipeline)
4. `SNV/Driver_Analysis/notebooks/OncoKB_Annotation_SNVs.ipynb` - filter SNVs to oncogenic & likely-oncogenic variants
5.  Run MutSig2CV (version 3.11) to identify driver genes (https://github.com/getzlab/MutSig2CV)
6. `SV/SV_sample_data/data_prep_for_oncoprint.ipynb` - prepare input file for SNV driver oncoprint
7. `SNV/Driver_Analysis/scripts/SNV_drivers_visualization.R` - SNV driver oncoprint (Extended Data Figure 1G)
