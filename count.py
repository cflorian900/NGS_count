'''
Perform operations to merge paired fastq files and count user defined sequences.
usage: python countSetup.py [--help, -h]
	positional arg (required): directory reference 
	Required args[[--names,-n][--user,-u]]
	Optional args[[--memory,-m][--processors,-p][--unzip,-gz][--merge,-m][--partitions,-parts][--names,-n][--user,-u]]
'''

import sys
import os
import pandas as pd
import numpy as np
import time
import dask.dataframe as dd
import concurrent.futures
import subprocess
from itertools import repeat
import processFQ
import argparse


def main(argv):
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('directory', help='Required: directory where data is stored (i.e. raw fastq files). Does not accept "." You must list directory tree.')
	parser.add_argument('reference', type=str, help='Required: file with two columns listing 1) sequence name and 2) sequence.')
	parser.add_argument('-mem','--memory', help='RAM per CPU per slurm job (defatult is 8).')
	parser.add_argument('-p','--processors', type=int, help='Request number of processors per slurm job (default is 1).')
	parser.add_argument('-gz','--unzip', help='Unzip paired fastq files: options are True or False (default is False).')
	parser.add_argument('-ta', '--trim3primeR1', type=str, help="Fasta file with 3' adapters to trim off of read 1.")
	parser.add_argument('-tg', '--trim5primeR1', type=str, help="Fasta file with 5' adapters to trim off of read 1.")
	parser.add_argument('-tA', '--trim3primeR2', type=str, help="Fasta file with 3' adapters to trim off of read 2.")
	parser.add_argument('-tG', '--trim5primeR2', type=str, help="Fasta file with 5' adapters to trim off of read 2.")
	parser.add_argument('-m','--merge', help='Use pandaseq to merge paired reads: options are True or False (default is True).')
	parser.add_argument('-parts','--partitions', type=int, help='Number of partitions to split up fastq file for counting (i.e. how many slurm jobs per fastq; default is 1).')
	parser.add_argument('-n','--names', type=str, help='File with two columns listing 1) sample name and 2) associated SIC indexes or some unique identifier in file name (no column headers).')
	parser.add_argument('-u','--user', type=str, help='Your HTCF username.')
	args = parser.parse_args()


	#for __init__ method
	#pFQ = processFQ.processFQ(args.directory)
	#fastqs = pFQ.generateFileList(args.directory)

	#for callable method
	pFQ = processFQ.processFQ()
	
	# generate the job array lookup file for cutadapt
	fastqs = pFQ.generateFileList(args.directory)

	# run cutadapt with user provided adapter fasta files
	# adapter trimming should be done before merging, but you need to fastqs object first
	# run cutadapt with user provided adapter fasta files
	# will probably need to make a trimmed directory, which will likely change some of the following code
	if any(item is not None for item in [args.trim3primeR1, args.trim5primeR1, args.trim3primeR2, args.trim5primeR2]):
		pFQ.trim(fastqs, args)
		# taking the sq approach described below to wait for job to finish
		sq = True
		while sq == True:
			out = subprocess.run('squeue -u %s' %(args.user),stdout=subprocess.PIPE, shell=True)
			# if len(out.stdout.decode('utf-8')) > 170: 
			if 'trim' in out.stdout.decode('utf-8'):
				time.sleep(10)	
			else:
				sq = False
	
	# change directory to trimmed but remember base directory
	baseDir = args.directory
	args.directory += '/trimmed/'

	# generate the job array lookup file for pandaseq
	trimFastqs = pFQ.generateFileList(args.directory)
	pFQ.lookupPandaseq(trimFastqs, baseDir, args)

	# run pandaseq job array; merged ouput is in the pandaseq directory
	pFQ.runPandaseq(args.directory, 'lookup_pandaseq.txt')

	# need to do something like this to wait for jobs to finish:
	# there's definitely a better way, but this is what we have for now
	sq = True
	while sq == True:
		out = subprocess.run('squeue -u %s' %(args.user),stdout=subprocess.PIPE, shell=True)
		# the stdout of squeue should be more than 85 characters if there is a job running
		# however, the driver job is always running 
		# an interactive job is ~158 characters, so I can safely say pandaseq is still running 
		# if the squeue character length is greater than 170 
		# if len(out.stdout.decode('utf-8')) > 170: 
		if 'panda' in out.stdout.decode('utf-8'):
			time.sleep(10)	
		else:
			sq = False

	# the parts below need to be in a different script submitted in the pandaseq.sh script to allow merging to finish
	# generate job array lookup file for counting merged data
	mergedFQ = pFQ.generateFileList(args.directory + 'pandaseq/')
	print('Files to process:', mergedFQ)
	numJobs = pFQ.countLookup(mergedFQ, args, baseDir)
	
	# write the analysis batch script 
	pFQ.countBatch(numJobs, args, baseDir)
	print('made counting batch script')

	# run analysis
	pFQ.runCount(args.directory)
	print('submitted counting batch script')

	sq = True
	while sq == True:
		out = subprocess.run('squeue -u %s' %(args.user),stdout=subprocess.PIPE, shell=True)
		# the stdout of squeue should be more than 85 characters if there is a job running
		# however, the driver job is always running 
		# an interactive job is ~158 characters, so I can safely say pandaseq is still running 
		# if the squeue character length is greater than 170 
		# if len(out.stdout.decode('utf-8')) > 170: #int(str(x.stdout.decode('utf-8'))[95:103]) > 0:
		if 'count_ar' in out.stdout.decode('utf-8') # count_ar is what shows up the squeue for running the counting portion
			time.sleep(60)	
		else:
			sq = False

	# make the merging batch script
	pFQ.makeMerge(args, baseDir)
	# merge the partitioned count files
	subprocess.run('sbatch %smerge.sh' %(args.directory + '/pandaseq/counts/'), shell=True)



if __name__ == '__main__':
	main(sys.argv)