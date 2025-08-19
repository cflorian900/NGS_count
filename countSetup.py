'''
Perform operations to merge paired fastq files and count user defined sequences.
usage: python countSetup.py [--help, -h]
	directory reference 
	Optional args: [--memory,-m][--processors,-p][--unzip,-gz][--merge,-m][--partitions,-parts][--names,-n][--user,-u]
'''

import sys
import argparse
import subprocess


def args(argv):
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('directory', help='Required: directory where data is stored (i.e. raw fastq files). Does not accept "." You must list directory tree.')
	parser.add_argument('reference', type=str, help='Required: file with two columns listing 1) sequence name and 2) sequence.')
	parser.add_argument('-mem','--memory', help='RAM per slurm job (defatult is 8).')
	parser.add_argument('-p','--processors', type=int, help='Request number of processors per slurm job (default is 1).')
	parser.add_argument('-gz','--unzip', help='Unzip paired fastq files: options are True or False (default is False).')
	parser.add_argument('-m','--merge', help='Use pandaseq to merge paired reads: options are True or False (default is True).')
	parser.add_argument('-parts','--partitions', type=int, help='Number of partitions to split up fastq file for counting (i.e. how many slurm jobs per fastq; default is 1).')
	parser.add_argument('-n','--names', type=str, help='File with two columns listing 1) sample name and 2) associated SIC indexes (no column headers).')
	parser.add_argument('-u','--user', type=str, help='Your HTCF username.')
	args = parser.parse_args()
	#print(args)

	with open('run_counter.sh', 'w') as f:
		f.write('#!/bin/bash\n#SBATCH --mem=1G\n#SBATCH --cpus-per-task=1\n' + 
			'#SBATCH -e count.err\n#SBATCH -o count.out\n\n' + 
			'cd %s\n\n' %(args.directory) +  
			'python count.py %s %s ' %(args.directory, args.reference) + 
			'-mem %s ' %(args.memory if args.memory else str(8)) + 
			'-p %s ' %(args.processors if args.processors else 1) + 
			'-gz %s ' %(args.unzip if args.unzip else False) + 
			'-m %s ' %(args.merge if args.merge else True) + 
			'-parts %s ' %(args.partitions if args.partitions else str(1)) + 
			'-n %s ' %(args.names if args.names else None) +
			'-u %s ' %(args.user))


args(sys.argv)

#if __name__ == '__main__':
#	args = args(sys.argv)
	# run the sbatch script that was just generated
#	subprocess.run('sbatch run_counter.sh %s' %(args.directory), shell=True)