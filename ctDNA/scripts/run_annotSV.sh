#!/usr/bin/env bash
#
# Run AnnotSV on structural variant calls
#
# Usage:
#     cd AnnotSV
#     ./run_annotsv.sh
#
# Inputs:
#   - ../data/sv_calls_v2.bed
#
# Outputs:
#   - ../results/ (main result is a .tsv)
#

set -euo pipefail

TOOLS_DIR="../../SV/tools"
ANNOTSV="${TOOLS_DIR}/AnnotSV"

INPUT="../data/leopard_cohort/annotsv_input.bed"
OUTDIR="../outputs/leopard_cohort/"

mkdir -p ${OUTDIR}

${ANNOTSV}/bin/AnnotSV \
  -SVinputFile ${INPUT} \
  -outputDir ${OUTDIR} \
  -genomeBuild GRCh38 \
  -annotationsDir ${ANNOTSV}/share/AnnotSV \
  -svtBEDcol 6 \
  -samplesidBEDcol 4
