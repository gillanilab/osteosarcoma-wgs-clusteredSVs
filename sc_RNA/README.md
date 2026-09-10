Region Expression
=================

### Workflow
1. `region_expression/notebooks/scRNA_analysis.ipynb` — identify subpopulation of cells within sample w/ hotspot clustered SV event , run deg analysis, prepare input files for cytotrace (Figure 4C, 4D; Extended Data Figure 2B, 2G)


cNMF Gene Expression Program (GEP) Analysis
============================================

### Workflow
1. `cnmf_gep_analysis/notebooks/cnmf_myc.ipynb` & `cnmf_gep_analysis/notebooks/cnmf_tp53.ipynb`— identify GEPs in MYC & TP53 upstream samples (Figure 4C, 4D, 4E; Extended Data Figure 2D) 
2. `cnmf_gep_analysis/notebooks/program_1_similarity` — test similarity between programs across the two samples (Extended Data Figure 2E)


CytoTRACE
=========

- `cytotrace/notebooks/cytotrace_all_samples.ipynb` - run cytoTRACE on all samples 
- `cytotrace/notebooks/cytotrace_myc.ipynb` & `cytotrace/notebooks/cytotrace_tp53.ipynb` - run cytoTRACE on individual samples of interest, input files generated from: `region_expression/notebooks/scRNA_analysis.ipynb` (Extended Data Figure 2C, 2G)


pySCENIC
========

### Workflow
1. `pyscenic/notebooks/prep-for-scenic.ipynb` 
2.  Run pySCENIC (https://github.com/aertslab/pySCENIC) 
3. `pyscenic/notebooks/scenic_analysis.ipynb` - associate regulons with programs (Extended Data Figure 2F)