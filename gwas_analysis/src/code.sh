#!/bin/bash

#run with --priority ---destnation rDNA\ Variations:GWAS --instance-type mem2_ssd1_v2_x4
#dx run -ichunk=0 -ichromosome=11 -y gwas_analysis --priority low --destination rDNA\ Variations:GWAS --instance-type mem1_hdd1_v2_x16

mkdir -p /home/dnanexus/output_folder

set -e -x
apt-get update
apt-get install -y --fix-missing wget tar unzip

chmod u+x /usr/bin/dxfuse

mkdir /home/dnanexus/mnt

#echo "Checking if mount works and we can see the input files"
dxfuse /home/dnanexus/mnt "rDNA Variations"

wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20250129.zip
unzip plink2_linux_avx2_20250129.zip

plink="./plink2"
plinkdir="/home/dnanexus/mnt/rDNA Variations/unrelated_WB"
mnt_dir="/home/dnanexus/mnt/rDNA Variations"

pgenfile=$plinkdir/"chr${chromosome}.pgen"
psamfile=$plinkdir/"chr${chromosome}.psam"
pvarfile=$plinkdir/"chr${chromosome}.pvar"
famfile=$plinkdir/"chr${chromosome}.fam"
bimfile=$plinkdir/"chr${chromosome}.bim"
covarfile=$mnt_dir/"ukbb_20pc.txt"
#
#ls -l $covarfile
#ls -l $phenofile

start=$((chunk * 1))
end=$((start + 1- 1))

samples_file="/home/dnanexus/mnt/rDNA Variations/run_list.txt"
# Read the CSV file and process the required chunk

sed -n "${start},${end}p" "${samples_file}" | while read -r variant_id;
do
  echo "$variant_id"
  phenofile=$mnt_dir/variant_phenotypes/"${variant_id}.txt"
  outdir=$variant_id
  mkdir -p $outdir
  outfile=$outdir/"chr${chromosome}"
  output_file=$outdir"/chr${chromosome}.PHENOTYPE.glm.linear"
  renamed_output=$outdir"/chr${chromosome}.${variant_id}.glm.linear"
  done_file="chr${chromosome}_${variant_id}.done"
  touch $done_file
  if [ -e "${outdir}/${done_file}" ]; then
    echo  "${renamed_output} File exists."
  else
    echo "Running Plink2 for chr${chromosome} ${variant_id}"
    $plink --pgen "$pgenfile" --bim "$bimfile" --fam "$famfile" --pheno "$phenofile" --covar "$covarfile" --glm hide-covar --covar-variance-standardize --out "$outfile" --threads 16 --vif 999
    mv "$output_file" "$renamed_output"
    output_files=()
    file_id=$(dx upload "${renamed_output}" --brief)
    dx-jobutil-add-output output_files "$file_id" --class=array:file
    file_id=$(dx upload "${done_file}" --brief)
    dx-jobutil-add-output output_files "$file_id" --class=array:file
  fi
done

