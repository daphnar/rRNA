# RGA.py file contains an implementation for the Reference Gap Alignment (RGA) algorithm.

The (RGA) steps are as followed:
1) We classified sequences as either 18S or 28S followed by Needleman–Wunsch global sequence alignment of each sequence to one RNA45S5 reference (either 18S or 28S based on read classification). 
2) We created a reference sequence that aligns to all other sequences that we call a gap-aligned reference. This gap-aligned reference has the same sequence as the reference, but at each nucleotide position, we extended a gap at the size of the maximal gap found by the global sequence alignment to all sequences. Importantly, this gap-aligned reference allows straightforward comparison among all sequences without requiring computationally expensive all-by-all pairwise sequence alignments. 
3) We aligned all H7-hESC sequences to the gap-aligned reference using the previous global alignment with additional extended gaps at reference positions. 
4) Lastly, we extracted all variants at a given position across all aligned sequences. 

Here we assume as input a fasta file where the first sequence in the file is the reference sequence. The following sequences will be aligned to the first sequence in the file.
 
