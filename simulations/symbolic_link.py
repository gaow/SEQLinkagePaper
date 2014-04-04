import sys
import subprocess

folder=sys.argv[1]

for gene in ['GJB2','SLC26A4','MYH9','MYO7A']:
	command='ln -s data/{}/{}* {}.tsv'.format(folder,gene,gene)
	subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE)
