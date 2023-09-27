# RGA.py file contains an implementation for the Reference Gap Alignment (RGA) algorithm.

## Expected input and returned output
Here we assume as input a fasta file where the first sequence in the file is the reference sequence. The following sequences will be aligned to the first sequence in the file.
When running the RGA.py the five sequences found in mock_sequences.fa file will be aligned to the first sequence, the 28S reference, and saved as the output file RGA_mock_sequences.fa

## System requirements
The system requirements are pre-installed Needleman-Wunch executable from https://www.ebi.ac.uk/Tools/emboss/. 
We use the executable "needle" from EMBOSS-6.6.0.
The code uses biopython and pandas packages.

## RGA method
The (RGA) steps are as follows where steps 1-4 are illustrated below:
1) We align using Needleman–Wunsch global sequence alignment each sequence to the reference sequence which is the first sequence in the provided fasta file.
2) We create a reference sequence that aligns with all other sequences that we call a gap-aligned reference. This gap-aligned reference has the same sequence as the reference, but at each nucleotide position, we extend a gap at the size of the maximal gap found by the global sequence alignment to all sequences. Importantly, this gap-aligned reference allows straightforward comparison among all sequences without requiring computationally expensive all-by-all pairwise sequence alignments. 
3) We align all sequences to the gap-aligned reference using the previous global alignment with additional extended gaps at reference positions. 
4) Lastly, we extract all variants at a given position across all aligned sequences. 
![Alt RGA](https://github.com/daphnar/rRNA/blob/main/RGA.png)




