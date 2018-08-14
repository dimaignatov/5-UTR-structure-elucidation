# 5-UTR-structure-elucidation

The software to analyse the FUSE (5' UTR Structure Eluscidation) data in L. monocytogenes

1) Load FATSTQ files from NCBI Sequence Reads Archive: SRP156446

2) Align the reads to L. monocytogenes EGD-e genome ('NC_003210.1.fa') with Bowtie2 aligner, using --end-to-end --very-sensitive mode.

3) Name the resulting files: 'FUSE_K.sam', 'FUSE_invitro.sam', 'FUSE_total.sam', 'FUSE_37A.sam', 'FUSE_37B.sam', 'FUSE_26.sam', 'FUSE_37shift.sam', 'FUSE_hfq.sam', 'FUSE_lhrA.sam',  'FUSE_prfA.sam'.

4) Run the script 'dms.py' in the folder containing these SAM files, L. monocytogenes  EGD-e genome sequence file 'NC_003210.1.fa' and the file 'NC_003210.1.gff3' containing genomic features. This script calculates statistics of mismatches, creates profiles for visualization in a genomic browser and the file 'Experiment.txt' containing information about coverage and number of mismatches for each nucleotide in the 5'UTRs and non-coding RNAs.
   
5) Calculate DMS values and perform comparison of DMS values with the script 'Test for changes of DMS values'.
   
