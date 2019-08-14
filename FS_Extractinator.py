#!/usr/bin/env python3

import argparse,sys,os,time
parser = argparse.ArgumentParser(prog='FS_Extractinator',usage='./FS_Extractinator.py [-h] [-i INPUT_FASTA] [-o OUTPUT_FASTA] [-f FLANKIG_SIZE] [-S START_POSITION] [-E END_POSITION] [-l SEQ_LENGTH] [-g] \n python3 FS_Extractinator.py [-h] [-i INPUT_FASTA] [-o OUTPUT_FASTA] [-f FLANKIG_SIZE] [-S START_POSITION] [-E END_POSITION] [-l SEQ_LENGTH] [-g]',
	description='''Extracts regions with flanking bases upstream and downstream for a specific position, given a fasta file. Applicable to multiple sequences.
	''',
								 epilog="Please send your suggestions and comments to Aditya  < lankapalli@shh.mpg.de >;) ")
parser.add_argument("-i", "--input",dest="input_fasta", help="Input multifasta")
parser.add_argument("-o", "--output",dest="output", help="Output prefix")
parser.add_argument("-l", "--length", dest="seq_length", help="Number of bases to be printed per line in output fasta [80]", type=int,default=80)
parser.add_argument("-f", "--flankingize", dest="flanking_size", help="Number of bases flanking upstream of the Start and downstream of the End position. Default size is 1000", type=int,default=1000)
parser.add_argument("-S", "--Start", dest="seq_Start", help="Start position of the region",type=int)
parser.add_argument("-E", "--End", dest="seq_End", help="End position of the region",type=int)

args = parser.parse_args()


with open(args.input_fasta) as f:
	g=list(filter(None,f.read().strip().split('>')))


align={}
for i in g:
	align[i.split('\n')[0]]=''.join(i.split('\n')[1:])


o=open(args.output+'.fasta','w')
u=open(args.output+'upstream.fasta','w')
d=open(args.output+'downstream.fasta','w')
for i in align:
	flanking_region=align[i][((args.seq_Start-1)-args.flanking_size):(args.seq_End+args.flanking_size)]
	flanking_up=align[i][((args.seq_Start-1)-args.flanking_size):(args.seq_Start-1)]
	flanking_down=align[i][args.seq_End:(args.seq_End+args.flanking_size)]
	o.write('>'+str(i)+'_'+str(((args.seq_Start)-args.flanking_size))+'_'+str((args.seq_End+args.flanking_size))+'\n')
	u.write('>'+str(i)+'_'+str(((args.seq_Start)-args.flanking_size))+'_'+str((args.seq_End+args.flanking_size))+'\n')
	d.write('>'+str(i)+'_'+str(((args.seq_Start)-args.flanking_size))+'_'+str((args.seq_End+args.flanking_size))+'\n')
	j=0
	for j in range(0,len(flanking_region),args.seq_length):
		o.write(str(flanking_region[j:j+args.seq_length])+'\n')
	k=0
	for k in range(0,len(flanking_up),args.seq_length):
		u.write(str(flanking_up[k:k+args.seq_length])+'\n')
	l=0
	for l in range(0,len(flanking_down),args.seq_length):
		d.write(str(flanking_down[l:l+args.seq_length])+'\n')
o.close()
u.close()
d.close()


sys.exit("{0} bases flanking a genomic region from {1} to {2} successfully written to {3}".format(args.flanking_size,args.seq_Start,args.seq_End,args.output))

