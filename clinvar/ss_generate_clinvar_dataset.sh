#!/bin/bash

if [ "$1" = "" ]
then
    echo -e "No argument about assembly version. Please add GRCh37 or GRCh38. \n"
    echo -e "Usage: "
    echo -e "        $0 [GRCh37|GRCh38] \n"
    exit 1
fi

set -euo pipefail
assembly="$1"
output_fn="clinvar_${assembly}.germline.nocoflicted.bcf.gz"
VCF="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${assembly}/clinvar.vcf.gz"
MD5="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${assembly}/clinvar.vcf.gz.md5"
TBI="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${assembly}/clinvar.vcf.gz.tbi"

function download_files() {
  echo "Donwloading ..."
  wget \
    --timestamping \
    --append-output=./wget_clnvcf.log \
    --show-progress \
    ${VCF}
  sleep 2
  md5sum clinvar.vcf.gz > downloaded_vcf.md5 &
  wget \
    --timestamping \
    --append-output=./wget_clnvcf.log \
    --show-progress \
    ${MD5}
  sleep 2
  wget \
    --timestamping \
    --append-output=./wget_clnvcf.log \
    --show-progress \
    ${TBI}
}

function check_md5() {
  echo "Check md5 ..."
  # md5sum clinvar.vcf.gz > downloaded_vcf.md5
  if md5sum --check downloaded_vcf.md5 | grep -Ev 'OK$' >/dev/null; then
    echo "The downloaded file may be corrupted."
    exit 2
  else
    echo "md5sum check -> OK! "
  fi
}

function filter_sort_bcf() {
  echo "Filtering ..."
  bcftools view \
    --output-type u \
    --include \
      'INFO/ORIGIN=="1" 
        & (INFO/CLNSIG=="Pathogenic" 
            | INFO/CLNSIG=="Likely_pathogenic" 
            | INFO/CLNSIG=="Pathogenic/Likely_pathogenic"
            | INFO/CLNSIG=="Benign" 
            | INFO/CLNSIG=="Likely_benign" 
            | INFO/CLNSIG=="Benign/Likely_benign")
        & INFO/CLNREVSTAT!="no_assertion_criteria_provided" 
        & INFO/CLNREVSTAT!="no_classification_provided" 
        & INFO/CLNREVSTAT!="no_classification_for_the_individual_variant"' \
    clinvar.vcf.gz \
  | bcftools sort \
      --output-type b \
      --output ${output_fn} \
      --write-index
}


function moving_files() {
  echo "Moving ..."
  local datetime=$(date +"%Y%m%d-%I%M%S")
  output_dir="Filtered_BCF_${assembly}_${datetime}"
  mkdir -p ${output_dir}
  mv ${output_fn}* ./${output_dir}
  mv wget_clnvcf.log ./${output_dir}
  # Remove downloaded files and md5
  rm -rf clinvar.vcf.gz* *.md5
}

function main() {
  download_files &&
  check_md5
  filter_sort_bcf
  moving_files
  echo -e "Completed!"
}


main


