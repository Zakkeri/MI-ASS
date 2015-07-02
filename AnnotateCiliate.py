# Used libraries
import argparse
import os.path
import re
import errno
import datetime
import time
from io import StringIO
from pyfasta import Fasta
from operator import itemgetter
from settings import *
from EssentialFunctions import *

# Debugging
DEBUGGING = True
time1 = time.time()

# Define argument parser
parser = argparse.ArgumentParser()
parser.add_argument('-mic', '--mic', dest='MIC')
parser.add_argument('-mac', '--mac', dest='MAC')
parser.add_argument('-o', '--o', dest='OUT')
parser.add_argument('-reblast', '--rb', dest='RB', action='store_true')

# Get arguments
if DEBUGGING:
	args = parser.parse_args('-mic ../Assembly_Data/Tetrahymena/tet_therm_processed-_mic_nuc.fa -mac ../Assembly_Data/Tetrahymena/Test_File.fasta -o ../Output_Tetrohymena'.split())
	#args = parser.parse_args('-mic ../Assembly_Data/Trifallax/oxy_tri_-_mic_assembly.fasta -mac ../Assembly_Data/Trifallax/Test_File.fasta -o ../Output_Trifallax --rb'.split())
	#args = parser.parse_args('-mic ./oxy_tri_-_mic_assembly.fasta -mac ./oxy_tri_-_mac_assembly_(with_pacbio).fasta -o ./Output'.split())
else:
	args = parser.parse_args()

# Store arguments in corresponding variables
regExp = Options['Tel_Reg_Exp']
MICfile = args.MIC
MACfile = args.MAC
Output_dir = args.OUT
ReBlast = args.RB

# Check if files are real files
while not os.path.isfile(MICfile):
	print(MICfile, ' is not a file, please input a valid MIC file:')
	MICfile = input()

while not os.path.isfile(MACfile):
	print(MACfile, ' is not a file, please input a valid MAC file:')
	MACfile = input()
	
# Test whether a regular expression is a valid one
re_comp = None
is_Valid_Re = False
while not is_Valid_Re:
	try:
		re_comp = re.compile(regExp, re.IGNORECASE)
		is_Valid_Re = True
	except re.error:
		print('The regular expression is invalid, please type a valid regular expression:')
		regExp = input()

# Create output directory if needed
LogFile = None
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
	LogFile.flush()
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
	output = subprocess.check_output("blastn -version".split(" "))
except:
	print('BLAST is not installed on the computer. Please install BLAST to run the program')
	logComment('BLAST is not installed on the computer, program terminated')
	exit()

logComment("BLAST version: " + output.decode(sys.stdout.encoding).split('\n')[0])

# Create BLAST database (and directory), if it doesn't exist
safeCreateDirectory(Output_dir + '/blast')
if not os.path.exists(Output_dir + '/blast/mic.nsq'):
	logComment('Building BLAST database from' + MICfile)
	output = subprocess.check_output(("makeblastdb -in " + str(MICfile) + ' -out ' +  Output_dir + '/blast/mic -dbtype nucl -parse_seqids -hash_index').split(" "))
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

#if DEBUGGING:
#	print(list(mac_fasta.keys()))

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

# Create output directory for MIC annotation
safeCreateDirectory(Output_dir + '/MIC_Annotation')	

# Create output directory for masked contig seqeunces
safeCreateDirectory(Output_dir + '/Masked_Contigs')	

# Create database input directory and database load files
if Options['DatabaseUpdate']:
	safeCreateDirectory(Output_dir + '/Database_Input')	
	temp = open(Output_dir + '/Database_Input/hsp.tsv', 'w')
	temp.close()
	temp = open(Output_dir + '/Database_Input/mds.tsv', 'w')
	temp.close()
	temp = open(Output_dir + '/Database_Input/tel.tsv', 'w')
	temp.close()

# Start BLASTing and annotating
logComment('Annotating ' + str(macCount) + ' MAC contigs...')

for contig in sorted(mac_fasta):
	if DEBUGGING:
		print("Annotating: ", str(contig))

	# create file for masked telomeres
	maskTel_file = open(Output_dir + '/Masked_Contigs/' + str(contig) + '.fa', 'w')
	maskTel_file.write(">" + str(contig) + "\n")
	
	# Run regular expression and mask telomeres
	seq = str(mac_fasta[contig])
	telomeres = re_comp.finditer(seq)
	left_Tel = None
	right_Tel = None
	
	# List for storing telomeric sequences that are not real telomeres
	tel_seq = list()
	
	# If regular expression was identified, then process
	if telomeres:
		str_pos = 0
		
		# Get start and end coordinates of every telomeric sequence
		tel_positions = [(m.span()[0], m.span()[1]) for m in telomeres]
		
		# Check for the telomere on the left
		indL = 0
		for coord in tel_positions:
			if coord[1] - coord[0] >= Options['TelomereLength'] and coord[0] <= Options['TelomereEndLimit']:
				maskTel_file.write(seq[str_pos:coord[0]].upper())
				maskTel_file.write(seq[coord[0]:coord[1]].lower())
				left_Tel = [coord[0],coord[1]]
				str_pos = coord[1]
				indL += 1
				break
			if coord[0] > Options['TelomereEndLimit']:
				break
			indL += 1
					
		# Check for telomere on the right
		for i in range(len(tel_positions) - 1, indL-1, -1):
			coord = tel_positions[i]
			if coord[1] - coord[0] >= Options['TelomereLength'] and len(mac_fasta[contig]) - coord[1] <= Options['TelomereEndLimit']:
				right_Tel = [coord[0], coord[1]]
				break
			if len(mac_fasta[contig]) - coord[1] > Options['TelomereEndLimit']:
				break
			
		# Mask telomeric sequence in the middle of contig
		for coord in tel_positions[indL:]:
			# If sequence is too short then skip it
			if coord[1] - coord[0] < Options['TelomereLength'] :
				continue
			# Otherwise, mask it
			if(str_pos != coord[0]):
				# Output seqeunce in upper case until the begining of the telomeric seqeunce
				maskTel_file.write(seq[str_pos:coord[0]].upper())
			maskTel_file.write(seq[coord[0]:coord[1]].lower())
			# If this is not the right telomere, then add it to the telomeric position
			if not right_Tel or right_Tel[0] != coord[0]:
				tel_seq.append((coord[0], coord[1]))
			str_pos = coord[1]
		
		# If there are still some letters left to print, print them
		if str_pos != len(seq):
			maskTel_file.write(seq[str_pos:len(seq)].upper())
		
		# check if left telomeres can be extended
		if left_Tel and tel_seq:
			numMerged = 0
			for tel in tel_seq:
				# If two telomeric seqeunces are within tolerance error, merge them
				if tel[0] - left_Tel[1] > Options["TelomericErrorTolerance"]:
					break
				left_Tel[1] = tel[1]
				numMerged += 1
			# Remove merged telomeres from the list
			if numMerged > 0:
				del tel_seq[:numMerged]
					
		# Check if right telomeres can be extended
		if right_Tel and tel_seq:
			numMerged = 0
			for tel in reversed(tel_seq):
				# If two telomeric seqeunces are within tolerance error, merge them
				if right_Tel[0] - tel[1] > Options["TelomericErrorTolerance"]:
					break
				right_Tel[0] = tel[0]
				numMerged += 1
			# Remove merged telomeres from the list
			if numMerged > 0:
				del tel_seq[-numMerged:]
							
	# Close masked contig file
	maskTel_file.close()
	
	# Run Rough BLAST pass
	MIC_maps = readBLAST_file(Output_dir + "/hsp/rough/" + str(contig) + ".csv") if not ReBlast and os.path.isfile(Output_dir + "/hsp/rough/" + str(contig) + ".csv") else run_Rough_BLAST(Output_dir, str(contig))
		
	# Get MAC start and MAC end with respect to telomeres
	#print(left_Tel, " ", right_Tel)
	MAC_start = 1 if not left_Tel else left_Tel[1]
	MAC_end = len(mac_fasta[contig]) if not right_Tel else right_Tel[0]
	#Debuging message		
	#print(contig, " start: ", MAC_start, " end: ", MAC_end)
	
	# Build list of MDSs if there is a BLAST output
	MDS_List = list()
	if MIC_maps != "":
		# Sort rough BLAST output
		sortHSP_List(MIC_maps)
		
		# Annotate MDSs with rough BLAST hsps
		getMDS_Annotation(MDS_List, MIC_maps, MAC_start, MAC_end)	
	else:
		MIC_maps = list()
	
	# Check if MAC is fully covered, and run fine pass if it is needed
	if getGapsList(MDS_List, MAC_start, MAC_end):
		# Run fine BLAST pass
		fineBLAST = readBLAST_file(Output_dir + "/hsp/fine/" + str(contig) + ".csv") if not ReBlast and os.path.isfile(Output_dir + "/hsp/fine/" + str(contig) + ".csv") else run_Fine_BLAST(Output_dir, str(contig))
		
		# If there is a BLAST output, process it
		if fineBLAST != "":
			# Sort fine BLAST output 
			sortHSP_List(fineBLAST)
			
			# Improve current annotation with fine BLAST results
			getMDS_Annotation(MDS_List, fineBLAST, MAC_start, MAC_end)
		
			# Add fine BLAST hsp into MIC_maps list if fine hsp is not a subset of some rough hsp
			temp_split = list()
			indexMap = {}
			ind = 0
			
			# Split hsps according to MIC contig name
			for hsp in MIC_maps:
				if hsp[1] in indexMap:
					temp_split[indexMap[hsp[1]]].append(hsp)
				else:
					indexMap[hsp[1]] = ind
					ind += 1
					temp_split.append([hsp])
			
			# Add fine BLAST results into rough blast list
			for hsp in fineBLAST:
				if hsp[1] not in indexMap or not [x for x in temp_split[indexMap[hsp[1]]] if int(x[5]) <= int(hsp[5]) and int(x[6]) >= int(hsp[6])]:
					MIC_maps.append(hsp)
			
	# Check for gaps and add them to the MDS List
	addGaps(MDS_List, MAC_start, MAC_end)	
		
	# Output results and label MDSs
	MDS_file = open(Output_dir + '/Annotated_MDS/' + str(contig) + '.tsv', 'w')
	ind = 1
	MDS_List.sort(key = lambda x: x[0])
	
	for mds in MDS_List[:-1]:
		MDS_file.write(str(ind) + '\t' + str(mds[0]) + '\t' + str(mds[1]) + '\t' + str(mds[2]) + '\n')
		mds.append(ind)
		ind += 1
	MDS_file.write(str(ind) + '\t' + str(MDS_List[-1][0]) + '\t' + str(MDS_List[-1][1]) + '\t' + str(MDS_List[-1][2]))
	MDS_List[-1].append(ind)
	MDS_file.close()
	
	# Annotate hsp with current MDS list
	mapHSP_to_MDS(MIC_maps, MDS_List)
	
	# Output MIC annotation
	MIC_maps.sort(key=lambda x: (x[1], int(x[7])))
	MIC_file = open(Output_dir + '/MIC_Annotation/' + str(contig) + '.tsv', 'w')
	prevMIC = ""
	for hsp in MIC_maps:
		if hsp[1] != prevMIC:
			MIC_file.write(hsp[1] + '\tMDS\tMAC start\tMAC end\tMIC start\tMIC end\n')
			prevMIC = hsp[1]
		MIC_file.write('\t' + str(hsp[-1]) + '\t' + hsp[5] + '\t' + hsp[6] + '\t' + hsp[7] + '\t' + hsp[8] + '\n')
		
	MIC_file.close()

	# Update database input files
	if Options['DatabaseUpdate']:
		updateDatabaseInput(MDS_List, MIC_maps, left_Tel, right_Tel, len(mac_fasta[contig]), Output_dir, contig)
	
logComment("Annotation is finished! Total time spent: " + str(time.time() - time1) + " seconds")
# Close all files
LogFile.close()

#Debugging time spent
print("Total time spent: ", time.time() - time1, " seconds")










