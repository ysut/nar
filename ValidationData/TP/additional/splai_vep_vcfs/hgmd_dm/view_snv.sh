#!/bin/zsh

set -eu

chr="$1"

bcftools view \
  --include 'STRLEN(REF)==1 && STRLEN(ALT)==1' \
  -Ov \
  all_DM_chr${chr}.splai.vep.vcf \
  -o hgmd_dm.chr${chr}.splai.vep.snv.vcf 
