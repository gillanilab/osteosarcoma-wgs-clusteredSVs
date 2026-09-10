# Clustered structural variant hotspots enable oncogenic addiction and plasticity in osteosarcoma

This repository contains the code used to analyze whole genome sequencing (WGS), bulk RNA sequencing, and single-cell RNA sequencing (scRNA seq) data from osteosarcoma patient tumor samples, as presented in the associated manuscript. Analysis is organized by data modality. Each subdirectory contains its own README with the analysis workflow, expected inputs and outputs, and any subdirectory specific dependencies.

### Repository Contents

| Directory | Description |
|---|---|
| `SNV` | De novo mutational signature analysis and artifact filtering (SignatureAnalyzer), variant annotation (Genome Nexus, OncoKB), and driver gene identification (MutSig2CV). |
| `SV` | Structural variant (SV) discovery and characterization, including SV landscape visualization, gene annotation (AnnotSV), recurrent breakpoint analysis (fishHook), chromothripsis detection (ShatterSeek), clustered SV and hotspot identification (signature.tools.lib), ecDNA and copy number analysis (AmpliconArchitect), and clonal analysis (SVclone). |
| `bulk_RNA` | Differential gene expression analysis of hotspot event-positive versus event-negative tumors. |
| `ctDNA` | Post-processing, clustered SV detection, and gene annotation (AnnotSV) of SV calls from ctDNA cohorts. |
| `sc_RNA` | Single-cell transcriptional program discovery and characterization, differentiation state inference (CytoTRACE), and regulon analysis (pySCENIC). |

Each subdirectory README documents the ordered workflow, the specific scripts and notebooks involved, and which manuscript figure each step produces where applicable.

## Data Availability

Raw sequencing data generated at Dana Farber Cancer Institute, Boston Children's Hospital, and the Broad Institute, including patient tumor WGS with matched normal controls and bulk RNA seq, will be deposited in dbGaP upon publication.

Whole genome sequencing data for pediatric tumor samples used in this study were obtained from St. Jude Cloud, dbGaP (phs000468, phs000699), and the European Genome phenome Archive (EGAD00001005389, EGAD00001002125).

Patient tumor snRNA seq data are available through the Alex's Lemonade Stand Foundation Single Cell Pediatric Cancer Atlas (ALSF scPCA), accessions SCPCP000017 and SCPCP000023. A processed version of the bulk patient tumor RNA seq data is also available under SCPCP000017.

### Citation

A manuscript is currently in preparation, and will be posted to a preprint server shortly!

### Contact

For questions, please contact the corresponding author Riaz Gillani (riaz_gillani@dfci.harvard.edu).
