#!/bin/bash

#run with --priority ---destnation rDNA\ Variations:output_by_mount
#dx run -ichunk=1 -y rdna_var_from_cram --priority low --destination rDNA\ Variations:output_by_mount
# Main script to download references, prepare environment, and run single_run.sh in parallel

# Reference files and directories
reference="Homo_sapiens_assembly38.fasta"
ribo_bed="ribo.bed"
bowtie2_index="bowtie2_atlas_expand150_ES/bowtie2_atlas_expand150_ES"
script="script_nuc_variant_freq_from_shortreads.py"
mkdir -p /home/dnanexus/output_folder

set -e -x
apt-get update
apt-get install -y --fix-missing samtools bowtie2 python3 git wget tar

chmod u+x /usr/bin/dxfuse
#dx-mount-all-inputs
#ls -l in

#wget https://go.dev/dl/go1.21.1.linux-amd64.tar.gz
#tar -C /usr/local -xzf go1.21.1.linux-amd64.tar.gz
#rm go1.21.1.linux-amd64.tar.gz
#export PATH="/usr/local/go/bin:${PATH}"
#
#git clone https://github.com/dnanexus/dxfuse.git
#cd dxfuse
#go build -o dxfuse cli/main.go

mkdir /home/dnanexus/mnt
#ls -l /home/dnanexus/dxfuse

#echo "Checking if mount works and we can see the reference genome"
dxfuse /home/dnanexus/mnt "rDNA Variations"

reference="/home/dnanexus/mnt/rDNA Variations/hg38_ref/Homo_sapiens_assembly38.fasta"
#sample="1041355_24048_0_0.dragen.cram"
#cram_file_path="/home/dnanexus/mnt/rDNA Variations/Bulk/DRAGEN WGS/Whole genome CRAM files (DRAGEN) [500k release]/$input_folder"
#sample_mount="$cram_file_path/$sample"
#echo "Value of folder_path: '${cram_file_path}'"

ls -l "$reference"

#chmod u+x ./single_run.sh
log_file="command_logs.txt"
touch $log_file

cram_file_path="/home/dnanexus/mnt/rDNA Variations/Bulk/DRAGEN WGS/Whole genome CRAM files (DRAGEN) [500k release]/"
past_run_output_path="/home/dnanexus/mnt/rDNA Variations/output_by_mount/"

count=0
output_files=()


start=$((chunk * 100))
end=$((start + 100 - 1))

samples_file="/home/dnanexus/mnt/rDNA Variations/folder_sample_UKBB.csv"
# Read the CSV file and process the required chunk

tail -n +2 "${samples_file}" | sed -n "${start},${end}p" | while IFS=',' read -r sample folder
do
    if [ -e "${past_run_output_path}/${sample}.counts.csv" ]; then
      echo "${sample}.counts.csv File exists." | tee -a $log_file
    else
      count=$((count + 1))
      echo "********Starting to work on sample number: ${count}"
      cram_file="$cram_file_path/$folder/$sample"

       # Step 1: Samtools view command to create .dragen.ribo.cram
      echo "Running samtools view for ${sample}" | tee -a $log_file
      samtools view -C -o "${sample}.ribo" --reference "$reference" -M -L "$ribo_bed" -@ 20 "${cram_file}"

      # Step 2: Samtools fastq to generate .dragen.ribo.fq.gz
      echo "Running samtools fastq for ${sample}" | tee -a $log_file
      samtools fastq -c 9 -0 /dev/null --reference "$reference" "${sample}.ribo" -o "${sample}.ribo.fq.gz"

      # Step 3: Bowtie2 alignment to create .dragen.ribo.atlas_mapped.sam
      echo "Running bowtie2 for ${sample}" | tee -a $log_file
      bowtie2 -x "$bowtie2_index" -U "${sample}.ribo.fq.gz" --score-min 'C,0,-1' -S "${sample}.ribo.atlas_mapped.sam" 2> "${sample}.mapping_stats.txt"
      grep -E "reads;|aligned" "${sample}.mapping_stats.txt" | awk '{print $1}' > "${sample}.alignment_counts.txt"
      #While accurate, this takes a long time - counting reads
      #read_count=$(samtools idxstats -@ 20 "${cram_file}" | awk '{sum += $3 + $4} END {print sum}')
      #Instead we can get the file size as an estimate of the number of reads
      read_count=$(stat -c%s "${cram_file}")
      read rDNA < <(sed -n '1p' "${sample}.alignment_counts.txt")
      read unmapped < <(sed -n '2p' "${sample}.alignment_counts.txt")
      read mapped_once < <(sed -n '3p' "${sample}.alignment_counts.txt")
      read mapped_multiple < <(sed -n '4p' "${sample}.alignment_counts.txt")
      echo -e "Total,$read_count\nrDNA,$rDNA\nUnmapped,$unmapped\nMapped_once,$mapped_once\nMapped_multiple,$mapped_multiple" > "${sample}.counts.csv"
      file_id=$(dx upload "${sample}.counts.csv" --brief)
      dx-jobutil-add-output output_files "$file_id" --class=array:file
#      file_id=$(dx upload "${sample}.alignment_counts.txt" --brief)
#      dx-jobutil-add-output output_files "$file_id" --class=array:file
      # Step 4: Run the custom Python script to generate .rdna_variants.csv
      echo "Running custom Python script for ${sample}" | tee -a $log_file
      python3 "$script" "${sample}.ribo.atlas_mapped.sam" "${sample}.rdna_variants.csv"
      echo "Created ${sample}.rdna_variants.csv" | tee -a $log_file
      file_id=$(dx upload "${sample}.rdna_variants.csv" --brief)
      dx-jobutil-add-output output_files "$file_id" --class=array:file
    fi
done

echo "All jobs finished."
done_file="done_${chunk}.txt"
touch $done_file
file_id=$(dx upload "${done_file}" --brief)
dx-jobutil-add-output output_files "$file_id" --class=array:file

