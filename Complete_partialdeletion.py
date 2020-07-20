#!/usr/bin/env python3
## Created by Aditya and modifications suggested by Felix Key
###USAGE: python3 completedeletioncode.py <fasta_filename> <percentage_of_deletion> <output_filename> <0 for variant sites only; anyother number above 0 for including invariant sites> '####
'''
Script filters sites that contain N in any seq in multifasta
4th argument set to 0, will filter in addition for variable sites only (i.e. >1 non-N character!)
'''

import argparse,sys
from collections import Counter

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
    Filters MultiFasta for CompleteDeletion (i.e. removes sites w/ "N"), and for invariable sites (if specified)
                                ''',
                                epilog="Questions or comments to Adytia ;) ")
parser.add_argument("-f", dest="fasta_file", help="Input multifasta", type=argparse.FileType('rt'))
parser.add_argument("-o", "--out_fasta", dest="out", help="Output multifasta", type=argparse.FileType('w'),default=sys.stdout)
parser.add_argument("-l", "--seq_length", dest="seq_length", help="Number of bases to be printed per line in output fasta [100]", type=int,default=100)
parser.add_argument("-p", "--percentage", dest="percentage", help="Percentage for partial deletion, default [100] is complete deletion ", type=int,default=100)
parser.add_argument('-n', dest="allSites", help="Set flag to print non-varaible sites too [False]", action='store_true', default=False) # default= not necessary as implied by action
args = parser.parse_args()

nucbases=['A','C','G','T']
percentage=args.percentage
out = args.out


## read the file and splits it and make a dictionary####
fastaseq=list(filter(None,open(args.fasta_file).read().strip().split('>')))
Alignment={i.split('\n')[0]:''.join(i.split('\n')[1:]) for i in fastaseq}


### Calculate the length of Alignment ###
p=set([len(Alignment[str(i)]) for i in Alignment])
if len(p)!=1:
	sys.exit(' All the sequences in the alignment are of not same length')
else:
	p=p.pop()

###select for seqeunces within the percentage of deletion ####







def var_invar(s):
    '''Identifies if a site is Variable or Invariable and makes a count of the bases'''
	F=[]
	basecomp=None
	F.append(s)
	basecomp=list(map(lambda x:Counter(s).__getitem__(x)/len(s),['A','C','G','T']))
	if len(set(s))>1:
		F.append('var')
	else:
		F.append('invar')
	F+=basecomp
	F+=[1-sum(basecomp),sum(basecomp)*100]
	return(F)



## build alignment matrix (lol) and filter for N (i.e. complete deletion)
n=n1=m1=[]
A=C=G=T=O=[0,0]
for j in range(0,p):
	s=[]
	w1=None
	for i in Alignment:
		s.append(Alignment[str(i)][j])
	w1=var_invar(s)
	if w1[-1]>=percentage:
		if w1[1]=='var':
			n.append(w1[0])
			A[0]=A[0]+w1[2]
			C[0]=C[0]+w1[3]
			G[0]=G[0]+w1[4]
			T[0]=T[0]+w1[5]
			O[0]=O[0]+w1[6]
			m1.append(w1[0])
		elif w1[1]=='invar':
			n1.append(w1[0])
			A[1]=A[1]+w1[2]
			C[1]=C[1]+w1[3]
			G[1]=G[1]+w1[4]
			T[1]=T[1]+w1[5]
			O[1]=O[1]+w1[6]
			m1.append(w1[0])
	else:
		pass


def rebuild(m1):
	Alignment1={}
	a=0
	for k in Alignment:
		n2=[]
		for l in m1:
			n2.append(l[a])
		Alignment1[k]=[''.join(n2)]
		a=a+1
	return(Alignment1)


if args.allSites==False:
	Alignment1=rebuild(n)
	print('printing alignment with only variable sites')
else:
	Alignment1=rebuild(m1)
	print('printing alignment with both variable and invariable sites')


##Reporting the Sequences and proportion of Nucleotides
print("Total Number of Sites in Alignment before deletion :",p)
print('\n')
print("Total Number of Sites in Alignment after "+str(percentage)+" deletion",len(m1))
print('\n')
print("Total Number of Variable sites : ",len(n))
print("Composition of A :",(A[0]/len(n)))
print("Composition of T :",(T[0]/len(n)))
print("Composition of G :",(G[0]/len(n)))
print("Composition of C :",(C[0]/len(n)))
print("Composition of Others :",(O[0]/len(n)))
print('\n\n')
print("Total Number of Invariable sites : ",len(n1))
print("Number of A sites:",A[1])
print("Number of T sites:",T[1])
print("Number of G sites:",G[1])
print("Number of C sites:",C[1])
print("Number of Other sites:",O[1])

### Write to the file ###
seq_length = args.seq_length
o = args.out
for i in Alignment1:
    # o.write('>'+str(i)+'\n'+Alignment1[str(i)][0]+'\n')
    o.write('>'+str(i)+'\n')
    sequence = Alignment1[str(i)][0]
    while len(sequence) > 0:
        o.write(sequence[:seq_length]+'\n')
        sequence = sequence[seq_length:]

o.close()
print("Successfully completed writing the alignment to ",args.out.name)
        

### END #####






