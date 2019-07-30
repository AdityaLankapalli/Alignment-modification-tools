# Alignment-modification-tools
Creates and Modifies genome-wide alignments
# Complete_Partialdeletion.py
Scripts to create alignments after deletion of gaps
<br>
requires Python>=3.6
<br>
chmod 755 completedeletioncode.py
<br>
completedeletioncode.py -f <fasta_filename> -p <percentage_of_deletion> -o <output_filename> -l <seq_wrap_length> -n [include non-variable sites or variable sites only]

{Also can be used as python3 completedeletioncode.py [options]}
<br>
<br>
Arguments in detail

-f  Input multifasta file
-o  Output prefix for multifasta file
-l  Number of bases to be printed per line in output fasta [100]
-p  Percentage for partial deletion, default [100] is for complete deletion
-n  Set flag to print non-varaible sites too [False]

