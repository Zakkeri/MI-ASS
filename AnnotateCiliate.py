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
#parser.add_argument('-re', '--re', dest='RE')
parser.add_argument('-mac', '--mac', dest='MAC')
parser.add_argument('-o', '--o', dest='OUT')

# Get arguments
if DEBUGGING:
	#args = parser.parse_args('-mic Test_Files/mic.fasta -mac Test_Files/mac.fasta -o ../Output'.split())
	args = parser.parse_args('-mic ../Assembly_Data/oxy_tri_-_mic_assembly.fasta -mac ../Assembly_Data/Test_File.fasta -o ../Output'.split())
	#args = parser.parse_args('-mic ./oxy_tri_-_mic_assembly.fasta -mac ./oxy_tri_-_mac_assembly_(with_pacbio).fasta -o ./Output'.split())
else:
	args = parser.parse_args()
	
regExp = Options['Tel_Reg_Exp']
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
	#call('blastn -version')
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

# Start BLASTing and annotating
logComment('Annotating ' + str(macCount) + ' MAC contigs...')

for contig in sorted(mac_fasta):
	if DEBUGGING:
		print("Annotating: ", str(contig))

	# create file for masked telomeres
	maskTel_file = open(Output_dir + '/hsp/rough/masked_' + str(contig) + '.fa', 'w')
	maskTel_file.write(">" + str(contig) + "\n")
	
	# Run regular expression and mask telomeres
	seq = str(mac_fasta[contig])
	telomeres = re_comp.finditer(seq)
	left_Tel = None
	right_Tel = None
	
	# List for storing telomeric sequences that are actually not telomeres
	tel_seq = list()
	
	# If regular expression was identified, then process
	if telomeres:
		str_pos = 0
		tel_positions = [(m.span()[0], m.span()[1]) for m in telomeres]
		# Check for the telomeres on the left
		indL = 0
		for coord in tel_positions:
			if coord[1] - coord[0] >= 10 and coord[0] <=100:
				maskTel_file.write(seq[str_pos:coord[0]].upper())
				maskTel_file.write(seq[coord[0]:coord[1]].lower())
				left_Tel = (coord[0],coord[1])
				str_pos = coord[1]
				indL += 1
				break
			if coord[0] > 100:
				break
			indL += 1
					
		# Check for telomere on the right
		for i in range(len(tel_positions) - 1, indL-1, -1):
			coord = tel_positions[i]
			if coord[1] - coord[0] >=10 and len(mac_fasta[contig]) - coord[1] <= 100:
				right_Tel = (coord[0], coord[1])
				break
			if len(mac_fasta[contig]) - coord[1] > 100:
				break
			
		# Mask telomeric sequence in the middle of contig
		for coord in tel_positions[indL:]:
			# If sequence is too short, or it is a right telomere, then skip it
			if coord[1] - coord[0] < 10 :
				continue
			# Otherwise, mask it
			if(str_pos != coord[0]):
				maskTel_file.write(seq[str_pos:coord[0]].upper())
			maskTel_file.write(seq[coord[0]:coord[1]].lower())
			# If this is not the right telomere, then add it to the telomeric position
			if not right_Tel or right_Tel[0] != coord[0]:
				tel_seq.append((coord[0], coord[1]))
			str_pos = coord[1]
		
		# If there are still some letters left to print, print them
		if str_pos != len(seq):
			maskTel_file.write(seq[str_pos:len(seq)].upper())
	
	# Close masked contig file
	maskTel_file.close()
	
	# Run Rough BLAST pass
	roughBLAST = run_Rough_BLAST(Output_dir, str(contig))
		
	# Sort rough BLAST output by 1) Coverage, 2) Length, 3) Persent  identity match, 4) MIC, 5) the hsp start position in the MAC
	MIC_maps = sorted(roughBLAST, key=lambda x: (float(x[11]), int(x[3]), float(x[2]), x[1], -int(x[5])), reverse=True)
	
	# Get MAC start and MAC end with respect to telomeres
	#print(left_Tel, " ", right_Tel)
	MAC_start = 1 if not left_Tel else left_Tel[1]
	MAC_end = len(mac_fasta[contig]) if not right_Tel else right_Tel[0]
	#Debuging message		
	#print(contig, " start: ", MAC_start, " end: ", MAC_end)
	
	# Build list of MDSs
	MDS_List = get_Rough_MDS_List(MIC_maps, MAC_start, MAC_end)	
	
	# Check if MAC is fully covered, and run fine pass if it is needed
	if getGapsList(MDS_List, MAC_start, MAC_end):
		# Run fine BLAST pass
		fineBLAST = run_Fine_BLAST(Output_dir, str(contig))
		if fineBLAST != "":
			# Improve current annotation with fine BLAST results
			improveAnnotation(fineBLAST, MDS_List, MAC_start, MAC_end)
		
			# Add fine BLAST hsp into MIC_maps list if fine hsp is not a subset of some rough hsp
			for hsp in fineBLAST:
				if not [x for x in MIC_maps if hsp[1] == x[1] and int(x[5]) <= int(hsp[5]) and int(x[6]) >= int(hsp[6])]:
					MIC_maps.append(hsp)
			
	# Check for gaps and add them to the MDS List
	addGaps(MDS_List, MAC_start, MAC_end)	
		
	# Output results and label MDSs
	MDS_file = open(Output_dir + '/Annotated_MDS/' + str(contig) + '.tsv', 'w')
	ind = 1
	MDS_List = sorted(MDS_List, key = lambda x: x[0])
	
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
	
logComment("Annotation is finished! Total time spent: " + str(time.time() - time1) + " seconds")
# Close all files
LogFile.close()

#Debugging time spent
print("Total time spent: ", time.time() - time1, " seconds")










