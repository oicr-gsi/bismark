#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

ls | sed 's/.*\.//' | sort | uniq -c
qrsh -l h_vmem=8G -N calculate_bismark_metrics -cwd -now n "\
. /oicr/local/Modules/default/init/bash; \
module load samtools 2>/dev/null; \
find . -regex '.*\.bam$' \
       -exec sh -c \" samtools flagstat {} | tr '\n' '\t'; echo \" \; \
| sort | uniq | tr '\t' '\n'"