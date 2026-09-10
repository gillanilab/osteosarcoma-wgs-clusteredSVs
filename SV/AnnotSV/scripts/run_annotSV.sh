#!/usr/bin/env bash
#
# Run AnnotSV on structural variant calls
#
# Usage:
#   From the repository root:
#     cd AnnotSV/scripts
#     ./run_annotsv.sh
#
# Inputs:
#   - ../data/sv_calls_v2.bed
#
# Outputs:
#   - ../results/ (main result is a .tsv)
#

set -euo pipefail

TOOLS_DIR="../../tools"
ANNOTSV="${TOOLS_DIR}/AnnotSV"

INPUT="../data/sv_calls_v2.bed"
OUTDIR="../results"

mkdir -p ${OUTDIR}

${ANNOTSV}/bin/AnnotSV \
  -SVinputFile ${INPUT} \
  -outputDir ${OUTDIR} \
  -genomeBuild GRCh38 \
  -annotationsDir ${ANNOTSV}/share/AnnotSV \
  -svtBEDcol 6 \
  -samplesidBEDcol 4
