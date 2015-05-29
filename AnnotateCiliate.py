# Used libraries
import argparse
import os.path
import re
import errno
import datetime
import subprocess
import sys
from io import StringIO
from subprocess import call
from pyfasta import Fasta
from operator import itemgetter
from settings import *
from EssentialFunctions import *

# Debugging
DEBUGGING = True

# Define argument parser
parser = argparse.ArgumentParser()
parser.add_argument('-mic', '--mic', dest='MIC')
parser.add_argument('-re', '--re', dest='RE')
parser.add_argument('-mac', '--mac', dest='MAC')
parser.add_argument('-o', '--o', dest='OUT')

# Get arguments
if DEBUGGING:
	args = parser.parse_args('-re AAAACCCC -mic Test_Files/mic.fasta -mac Test_Files/mac.fasta -o ../Output'.split())
else:
	args = parser.parse_args()
	
regExp = args.RE
MICfile = args.MIC
MACfile = args.MAC
Output_dir = args.OUT

# Check if files are real files
while not os.path.isfile(MICfile):
	print(MICfile, ' is not a file, please input a valid MIC file:')
	MICfile = input()

while not os.path.isfile(MACfile):
	print(MACfile, ' is not a file, please input a valid MAC file:')
	MACfile = input()
	
# Test whether a regulat expression is a valid one
re_comp = 0
is_Valid_Re = False
while not is_Valid_Re:
	try:
		re_comp = re.compile(regExp)
		is_Valid_Re = True
	except re.error:
		print('The regular expression is invalid, please type a valid regular expression:')
		regExp = input()

LogFile = None
# Create output directory if needed
try:
	os.makedirs(Output_dir)
	LogFile = open(Output_dir + '/log.txt', 'w')
	LogFile.write(datetime.datetime.now().strftime("%I:%M%p %B %d %Y") + ' - Directory ' + Output_dir + ' created\n')
except OSError as exception:
	if exception.errno == errno.EEXIST:
		LogFile = open(Output_dir + '/log.txt', 'w')
		LogFile.write(datetime.datetime.now().strftime("%I:%M%p %B %d %Y") + ' - Directory ' + Output_dir + ' found\n')
	else:
		raise

#-----------------------------------------------------------------------------------------------------
# Define log comment function for the ease of use
def logComment(comment):
	global LogFile
	LogFile.write(datetime.datetime.now().strftime("%I:%M%p %B %d %Y") + ' - ' + comment + '\n')
#-----------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------
# Define create directory function
def safeCreateDirectory(dir):
	try:
		os.makedirs(dir)
		logComment('Directory ' + dir + ' created')
	except OSError as exception:
		if exception.errno != errno.EEXIST:
			raise
#-----------------------------------------------------------------------------------------------------
		
# Print parsed arguments
if DEBUGGING:
	print(regExp, MICfile, MACfile)

# Check if BLAST is installed
output = None
try:
	#call('blastn -version')
	output = subprocess.check_output("blastn -version")
except:
	print('BLAST is not installed on the computer. Please install BLAST to run the program')
	logComment('BLAST is not installed on the computer, program terminated')
	exit()

logComment("BLAST version: " + output.decode(sys.stdout.encoding).split('\n')[0])

# Create BLAST database (and directory), if it doesn't exist
safeCreateDirectory(Output_dir + '/blast')
if not os.path.exists(Output_dir + '/blast/mic.nsq'):
	logComment('Building BLAST database from' + MICfile)
	output = subprocess.check_output("makeblastdb -in " + str(MICfile) + ' -out ' +  Output_dir + '/blast/mic -dbtype nucl -parse_seqids -hash_index')
	logComment(output.decode(sys.stdout.encoding))
else:
	logComment("BLAST database for " + MICfile + " found")

# Import MAC fasta file
logComment('Importing MAC fasta file...')
mac_fasta = None
try:
	mac_fasta = Fasta(MACfile)
except Exception as e:
	print("Error while importing fasta file\n" + str(e))
	logComment("Can't import fasta file\n" + str(e))
	exit()

if DEBUGGING:
	print(list(mac_fasta.keys()))

# Record number of imported MAC contigs
macCount = len(mac_fasta.keys())
logComment(str(macCount) + ' sequences imported')

# Create hsp, hsp/rough, hsp/fine directtories
safeCreateDirectory(Output_dir + '/hsp')
safeCreateDirectory(Output_dir + '/hsp/rough')
safeCreateDirectory(Output_dir + '/hsp/fine')

# Rough Blast parameters
dust = "yes" if Options['RoughBlastDust'] else "no"
ungapped = " -ungapped " if Options['RoughBlastUngapped'] else ""
maskLowercase = " -lcase_masking " if Options['BlastMaskLowercase'] else ""
	
logComment("BLAST rough pass parameters:\nblastn -task " + Options['RoughBlastTask'] + " -word_size " + str(Options['RoughBlastWordSize']) + " -max_hsps 0 " +
"-max_target_seqs 10000 -dust " + dust + ungapped + maskLowercase + "-num_threads " + str(Options['ThreadCount']) + 
" -outfmt \"10 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qcovs\"\n")

#Fine Blast parameters
dust = "yes" if Options['FineBlastDust'] else "no"
ungapped = " -ungapped " if Options['FineBlastUngapped'] else ""
	
logComment("BLAST fine pass parameters:\nblastn -task " + Options['FineBlastTask'] + " -word_size " + str(Options['FineBlastWordSize']) + " -max_hsps 0 " + 
"-max_target_seqs 10000 -dust " + dust + ungapped + maskLowercase + 
"-outfmt \"10 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qcovs\"" + " -num_threads " + str(Options['ThreadCount']) + "\n")

# Create output directory for MAC mds
safeCreateDirectory(Output_dir + '/Annotated_MDS')	

# Start BLASTing and annotating
logComment('Annotating ' + str(macCount) + ' MAC contigs...')

for contig in mac_fasta:
	# create file for masked telomeres
	maskTel_file = open(Output_dir + '/hsp/rough/masked_' + str(contig) + '.fa', 'w')
	maskTel_file.write(">" + str(contig) + "\n")
	
	# Run regular expression and mask telomeres
	seq = str(mac_fasta[contig])
	telomeres = re_comp.finditer(seq)
	tell_pos = list()
	str_pos = 0
	for iter in telomeres:
		coord = iter.span()
		if(str_pos != coord[0]):
			maskTel_file.write(seq[str_pos:coord[0]].upper())
		maskTel_file.write(seq[coord[0]:coord[1]].lower())
		tell_pos.append((coord[0], coord[1]))
		str_pos = coord[1]
	if str_pos != len(seq):
		maskTel_file.write(seq[str_pos:len(seq)].upper())
	
	# Close masked contig file
	maskTel_file.close()
	
	# Run Rough BLAST pass
	dust = "yes" if Options['RoughBlastDust'] else "no"
	ungapped = " -ungapped " if Options['RoughBlastUngapped'] else ""
	maskLowercase = " -lcase_masking " if Options['BlastMaskLowercase'] else ""

	rough_out = subprocess.check_output("blastn -task " + Options['RoughBlastTask'] + " -word_size " + str(Options['RoughBlastWordSize']) + " -max_hsps 0 " +
	"-max_target_seqs 10000 -dust " + dust + ungapped + maskLowercase + "-query " + Output_dir + "/hsp/rough/masked_" + str(contig) + ".fa -db " +
	Output_dir + "/blast/mic -num_threads " + str(Options['ThreadCount']) + 
	" -outfmt \"10 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qcovs\"")
	
	# Filter empty rows
	roughVal = [x.rstrip() for x in rough_out.decode(sys.stdout.encoding).split('\n') if x != ""]
		
	# Save rough results to file
	rough_file = open(Output_dir + '/hsp/rough/' + str(contig) + '.csv', 'w')
	for x in roughVal[:-1]:
		rough_file.write(x + '\n')
	rough_file.write(roughVal[-1])	
	rough_file.close()
	
	# Run Fine BLAST pass
	dust = "yes" if Options['FineBlastDust'] else "no"
	ungapped = " -ungapped " if Options['FineBlastUngapped'] else ""

	fine_out = subprocess.check_output("blastn -task " + Options['FineBlastTask'] + " -word_size " + str(Options['FineBlastWordSize']) + " -max_hsps 0 " + 
	"-max_target_seqs 10000 -dust " + dust + ungapped + maskLowercase + "-query " + Output_dir + "/hsp/rough/masked_" + str(contig) + ".fa " +
	"-db " + Output_dir + "/blast/mic " +
	"-outfmt \"10 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qcovs\"" + " -num_threads " + str(Options['ThreadCount']))
	
	# Filter empty rows
	fineVal = [x.rstrip() for x in fine_out.decode(sys.stdout.encoding).split('\n') if x != ""]
	
	# Save fine results to file
	fine_file = open(Output_dir + '/hsp/fine/' + str(contig) + '.csv', 'w')
	for x in fineVal[:-1]:
		fine_file.write(x + "\n")
	fine_file.write(fineVal[-1])
	fine_file.close()
	
	# Combine result of rough and fine blast and remove duplicates
	res = [x.split(',') for x in list(set(roughVal).union(set(fineVal)))]
		
	# Split list according to MIC contigs and sort every splitted entry by the hsp start position in the MAC
	MIC_maps = list()
	for mic in {x[1] for x in res}:
		temp = sorted([hsp for hsp in res if hsp[1] == mic], key=lambda x: int(x[5]))
		MIC_maps.append(temp)
	
	# Get MAC start and MAC end with respect to telomeres
	MAC_start = 0
	MAC_end = len(mac_fasta[contig])
	if len(tell_pos) >= 2:
		MAC_start = int(tell_pos[0][1])
		MAC_end = int(tell_pos[-1][0])
	elif len(tell_pos) == 1:
		if int(tell_pos[0][0]) < MAC_end - int(tell_pos[0][1]):
			MAC_start = int(tell_pos[0][1])
		else:
			MAC_end = int(tell_pos[0][0])
	#Debuging message		
	#print(contig, " start: ", MAC_start, " end: ", MAC_end)
	
	# Build list of MDSs
	MDS_List = getMDS_List(MIC_maps)	
	
	# Check for gaps and add them to the MDS List
	addGaps(MDS_List, MAC_start, MAC_end)	
		
	# Output results
	MDS_file = open(Output_dir + '/Annotated_MDS/' + str(contig) + '.tsv', 'w')
	ind = 1
	MDS_List = sorted(MDS_List, key = lambda x: x[0])
	
	for mds in MDS_List[:-1]:
		MDS_file.write(str(ind) + '\t' + str(mds[0]) + '\t' + str(mds[1]) + '\t' + str(mds[2]) + '\n')
		ind += 1
	MDS_file.write(str(ind) + '\t' + str(MDS_List[-1][0]) + '\t' + str(MDS_List[-1][1]) + '\t' + str(MDS_List[-1][2]))
	MDS_file.close()
	
	
# Close all files
LogFile.close()












