# RiboVAn.py file contains an implementation for the RiboVAn pipeline.

The RiboVAn.py script expect as input parameters mapped sam files to the atlas. 

Here since input files were paired-end short reads, both strand reads were mapped to the atlas. 
As an example:
fastq_file1="file_1.fastq.gz"
fastq_file2="file_2.fastq.gz"

bowtie2 -x bowtie2_inde} -U "$fastq_file1" --score-min 'C,0,-1' \
-S "file_1.ribo.atlas_mapped.sam" --un "file_1.unmapped.sam" \
2> "file_1.mapping_stats.txt"

bowtie2 -x bowtie2_inde} -U "$fastq_file1" --score-min 'C,0,-1' \
-S "file_2.ribo.atlas_mapped.sam" --un "file_1.unmapped.sam" \
2> "file_2.mapping_stats.txt"

Then file_2.mapping_stats.txt and file_2.mapping_stats.txt are entered as input files to RiboVAn.py as followed:

python RiboVAn.py file_1.ribo.atlas_mapped.sam file_2.ribo.atlas_mapped.sam output

output result will contain nucleodite atlas variant frequencies.
