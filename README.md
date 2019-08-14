# Alignment-modification-tools
_Create and Modify genome-wide alignments_

# Complete_Partialdeletion.py
Scripts to create alignments after deletion of gaps
<br>
### requirements
requires Python>=3.6
<br>
### change permissions
chmod 755 completedeletioncode.py
<br>
```
completedeletioncode.py -f <fasta_filename> -p <percentage_of_deletion> -o <output_filename> -l <seq_wrap_length> -n [include non-variable sites or variable sites only]
```
### _Also can be used as python3 completedeletioncode.py [options]_
<br>
<br>
**Arguments in detail**

-f  Input multifasta file<br>
-o  Output prefix for multifasta file<br>
-l  Number of bases to be printed per line in output fasta [100]<br>
-p  Percentage for partial deletion, default [100] is for complete deletion<br>
-n  Set flag to print non-varaible sites too [False]<br>


# FS_Extractinator.py
Script to extract upstream and downstream sequences for a specificed region
<br>
### requirements
requires Python>=3.6
<br>
### change permissions
chmod 755 FS_Extractinator.py
<br>
```
FS_Extractinator.py -i <multi_fasta_filename> -o <output_prefix> -l <seq_wrap_length> -S <Start position> -E <End position> -f <size of up/down flanking the sequence positions>
 ```

**Arguments in detail**

-i  Input multifasta file<br>
-o  Output prefix<br>
-l  Number of bases to be printed per line in output fasta [80]<br>
-f  Number of bases upstream and downstream flanking the Start and the End position. Default size is 1000<br>
-S  Start position of the region<br>
-E  End position of the region<br>
