# Python file with all essential functions used during CiliateAnnotation program execution
import subprocess
import sys
import datetime
import time
import os.path
import errno
import math
import regex as re
from datetime import datetime
from settings import *
from functools import reduce

#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Define log comment function for the ease of use

def logComment(comment):
	logComment.logFile.write(datetime.now().strftime("%I:%M%p %B %d %Y") + ' - ' + comment + '\n')
	logComment.logFile.flush()

# Log file
logComment.logFile = None
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Define create directory function

def safeCreateDirectory(dir):
	try:
		os.makedirs(dir)
		logComment('Directory ' + dir + ' created')
	except OSError as exception:
		if exception.errno != errno.EEXIST:
			raise
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function creates all necessary output folders and files

def createOutputDirectories(Output_dir):
	# Create BLAST database directory
	safeCreateDirectory(Output_dir + '/blast')

	# Create hsp, hsp/rough, hsp/fine directtories
	safeCreateDirectory(Output_dir + '/hsp')
	safeCreateDirectory(Output_dir + '/hsp/rough')
	safeCreateDirectory(Output_dir + '/hsp/fine')
	
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
		temp = open(Output_dir + '/Database_Input/arr.tsv', 'w')
		temp.close()
	
	# Create output gff3 directory and file
	safeCreateDirectory(Output_dir + '/GFF')
	gffFile = open(Output_dir + "/GFF/mac_mds.gff3", "w")
	gffFile.write("##gff-version 3\n")
	gffFile.close()
	
	# Create output directory and files for MIC scrambling patterns
	safeCreateDirectory(Output_dir + "/Scrambling")
	temp = open(Output_dir + "/Scrambling/all.tsv", "w")
	temp.close()
	temp = open(Output_dir + "/Scrambling/scrambled.tsv", "w")
	temp.close()
	temp = open(Output_dir + "/Scrambling/maps.tsv", "w")
	temp.close()
	
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function runs rough BLAST and returns the hsp result list

def run_Rough_BLAST(Output_dir, contig):
	# Set up parameters
	param = ("blastn -task " + Options['RoughBlastTask'] + " -word_size " + str(Options['RoughBlastWordSize']) + " -max_hsps 0 -max_target_seqs 10000 -dust").split(" ")
	
	if Options['RoughBlastDust']:
		param.append("yes")
	else:
		param.append("no")
	if Options['RoughBlastUngapped']:
		param.append("-ungapped")
	if Options['BlastMaskLowercase']:
		param.append("-lcase_masking")

	param.append("-query")
	param.append(Output_dir + "/Masked_Contigs/" + str(contig) + ".fa")
	param.append("-db")
	param.append(Output_dir + "/blast/mic")
	param.append("-num_threads")
	param.append(str(Options['ThreadCount']))
	param.append("-outfmt")
	param.append("10 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qcovs")	

	# Run BLAST command
	rough_out = subprocess.check_output(param)
	
	# Filter empty rows
	roughVal = [x.rstrip() for x in rough_out.decode(sys.stdout.encoding).split('\n') if x != "" and float(x.rstrip().split(",")[11]) >= Options['MIC_Coverage_Threshold']]
	# Check if there are any hsps
	if not roughVal:
		return ""
	
	# Save rough results to file
	rough_file = open(Output_dir + '/hsp/rough/' + contig + '.csv', 'w')
	for x in roughVal[:-1]:
		rough_file.write(x + '\n')
	rough_file.write(roughVal[-1])	
	rough_file.close()
	
	# Filter duplicates out and parse result
	res = [x.split(',') for x in list(set(roughVal))]
	
	return res
	
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function runs fine BLAST and returns hsp result list

def run_Fine_BLAST(Output_dir, contig):
	param = ("blastn -task " + Options['FineBlastTask'] + " -word_size " + str(Options['FineBlastWordSize']) + " -max_hsps 0 -max_target_seqs 10000 -dust").split(" ")
 	
	if Options['FineBlastDust']:
		param.append("yes")
	else:
		param.append("no")
	if Options['FineBlastUngapped']:
		param.append("-ungapped")
	if Options['BlastMaskLowercase']:
		param.append("-lcase_masking")
	
	param.append("-query")
	param.append(Output_dir + "/Masked_Contigs/" + str(contig) + ".fa")
	param.append("-db")
	param.append(Output_dir + "/blast/mic")
	param.append("-num_threads")
	param.append(str(Options['ThreadCount']))
	param.append("-outfmt")
	param.append("10 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qcovs")

	fine_out = subprocess.check_output(param)
	
	# Filter empty rows
	fineVal = [x.rstrip() for x in fine_out.decode(sys.stdout.encoding).split('\n') if x != "" and float(x.rstrip().split(",")[11]) >= Options['MIC_Coverage_Threshold']]
	# Check if there are any hsps
	if not fineVal:
		return ""
	
	# Save fine results to file
	fine_file = open(Output_dir + '/hsp/fine/' + str(contig) + '.csv', 'w')
	for x in fineVal[:-1]:
		fine_file.write(x + "\n")
	fine_file.write(fineVal[-1])
	fine_file.close()
	
	# Filter duplicates out and parse result
	res = [x.split(',') for x in list(set(fineVal))]
	
	return res
	
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function reads hsp file that was generated by run_Rough_BLAST or run_Fine_BLAST functions previously
def readBLAST_file(Filename):
	# Initialize list to read input
	res = list()
	
	# Read file line by line
	for line in open(Filename, "r"):
		# Parse and append line into the res list
		parsed =  line.rstrip().split(",")
		if float(parsed[11]) >= Options['MIC_Coverage_Threshold']:
			res.append(parsed)
		
	# Return list if it not empty and empty string otherwise
	if res:
		return res
	else:
		return ""
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function goes through the list of high scoring pairs that are associated with the MIC and constructs MDSs for the MAC

def getMDS_Annotation(MDS_List, HSP_List, MAC_start, MAC_end):
	# Get Gaps List
	Gaps = getGapsList(MDS_List, MAC_start, MAC_end)
		
	# Build list of MDSs
	for hsp in HSP_List:
		mds_toAdd = [int(hsp[5]), int(hsp[6]),0]
				
		# If current MDS does not overlap with any gap, then skip it
		if not [x for x in Gaps if x[1] > mds_toAdd[0] and x[0] < mds_toAdd[1]]:
			continue
		
		# Go through MDSs that overlap with current MDS and see if any can be merged or removed
		for x in sorted([mds for mds in MDS_List if mds[1] > mds_toAdd[0] and mds[0] < mds_toAdd[1]]):
			# If x is a subset of mds_toAdd, then remove x
			if x[0] >= mds_toAdd[0] and x[1] <= mds_toAdd[1]:
				MDS_List.remove(x)
			# Check if two MDSs can be merged
			elif int((mds_toAdd[0] + mds_toAdd[1])/2) in range(x[0], x[1]) or int((x[0] + x[1])/2) in range(mds_toAdd[0], mds_toAdd[1]):
				mds_toAdd[0] = min(mds_toAdd[0], x[0])
				mds_toAdd[1] = max(mds_toAdd[1], x[1])
				MDS_List.remove(x)
		
		# Add MDS to the MDS list, update gaps list, check if we are done
		MDS_List.append(mds_toAdd)
		Gaps = getGapsList(MDS_List, MAC_start, MAC_end)
		if not Gaps:
			break
	
	# Sort the MDS List
	MDS_List.sort(key=lambda x: x[0])

	
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This function takes a list of MDSs (sorted by the MDS begining coordinate) and returns the intervals of the MAC covering

def getCovering_Intervals(MDS_List):
	# Construct intervals by using MDS List and checking for gaps between consecutive MDSs
	MAC_Interval = list()
	for mds in sorted(MDS_List, key = lambda x: x[0]):
		if not MAC_Interval:
			MAC_Interval.append([mds[0], mds[1]])
		else:
			if MAC_Interval[-1][1] >= mds[0] - 1:
				MAC_Interval[-1][1] = max(mds[1], MAC_Interval[-1][1])
			else:
				MAC_Interval.append([mds[0], mds[1]])
	
	return MAC_Interval
				
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function checks whether there are any gaps in the MAC annotation

def addGaps(MDS_List, MAC_start, MAC_end):
	# If MDS List is empty, then return the whole MAC interval as a gap
	if not MDS_List:
		MDS_List.append([MAC_start, MAC_end, 1])
		return
	
	# Get the list of covered MAC Interval(s)
	MAC_Interval = getCovering_Intervals(MDS_List)
		
	# If we have more than one interval, then there are gaps we need to add to the annotation
	if len(MAC_Interval) > 1:
		prev = MAC_Interval[0]
		for interv in MAC_Interval[1:]:
			MDS_List.append([prev[1], interv[0],1])
			prev = interv
			
	# Check for gaps at the begining of MAC and at the end of MAC
	if MAC_Interval[0][0] - MAC_start > 0:
		MDS_List.append([MAC_start, MAC_Interval[0][0], 1])
	if MAC_end - MAC_Interval[-1][1] > 0:
		MDS_List.append([MAC_Interval[-1][1], MAC_end, 1])
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function assigns MDS number to hsps that correspond to some MDS in MAC

def mapHSP_to_MDS(MIC_maps, MDS_List):
	# Go through the list of hsps and assign MDS number, or -1 to each hsp
	for hsp in MIC_maps:
		# Set hsp to no MDS for now
		hsp.append(-1)
			
		# Get list of MDSs that were mapped from current hsp
		overlap = [x for x in MDS_List if (x[2] != 1) and (int(hsp[5]) < x[1] and int(hsp[6]) > x[0])]
		if not overlap:
			continue
			
		# Define reduce function to decide what MDS the hsp is going to match the best (has biggest overlap)
		match = lambda a, b: a if min(a[1], int(hsp[6])) - max(a[0], int(hsp[5])) > min(b[1], int(hsp[6])) - max(b[0], int(hsp[5])) else b
		matched_MDS = reduce(match, overlap)
		
		# Assign hsp to MDS
		hsp[-1] = matched_MDS[-1]
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function calculates the list of gaps in the MAC annotation
		
def getGapsList(MDS_List, MAC_start, MAC_end):
	# If no MDS, return the whole contig interval
	if not MDS_List:
		return [[MAC_start, MAC_end]]
	
	# Sort MDS list
	MDS_List.sort(key=lambda x: x[0])
	
	# Add gap at the begining, if needed
	Gaps = list()
	if MDS_List[0][0] - MAC_start > 1:
		Gaps.append([MAC_start, MDS_List[0][0]])
	
	# If there are more than one MDS, then add gaps in between MDSs
	if len(MDS_List) > 1:
		prev = MDS_List[0]
		for x in MDS_List[1:]:
			if x[0] - prev[1] > 1:
				Gaps.append([prev[1], x[0]])
			prev = x
	
	# Check if there is a gap at the end	
	if MAC_end - MDS_List[-1][1] > 1:
		Gaps.append([MDS_List[-1][1], MAC_end])
		
	return Gaps
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function outputs annotation results into the database load file

def updateDatabaseInput(MDS_List, MIC_maps, MIC_to_HSP, left_Tel, right_Tel, mac_length, Output_dir, contig):
	# If MIC_maps are not empty, then update hsp file
	if MIC_maps:
		# Open hsp file to append
		hspFile = open(Output_dir + '/Database_Input/hsp.tsv', 'a')
		
		# Output hsps to file
		for hsp in MIC_maps:
									
			# Get MIC start, end, and orientation
			micStart = hsp[7]
			micEnd = hsp[8]
			micOrient = "+"
			if int(hsp[7]) > int(hsp[8]):
				micStart = hsp[8]
				micEnd = hsp[7]
				micOrient = "-"
			
			# Print hsp
			hspFile.write("\\N\t" + hsp[0] + "\t\\N\t" + str(hsp[-1]) + "\t" + hsp[5] + "\t" + hsp[6] + "\t" + hsp[1] + 
			"\t\\N\t" + micStart + "\t" + micEnd + "\t" + micOrient + "\t" + hsp[3] + "\t" + hsp[2] + "\t" + hsp[4] + "\t" + hsp[9] + "\t" + hsp[10] + 
			"\t" + hsp[11] + "\n")
		
		# close hsp file
		hspFile.close()
	
	# If MDS list is not empty, then update mds file
	if MDS_List:
		# Open mds file to append
		mdsFile = open(Output_dir + '/Database_Input/mds.tsv', 'a')
		
		# Output mdss to file
		for mds in MDS_List:
			#Print mds
			mdsFile.write("\\N\t\\N\t" + contig + "\t" + str(mds[-1]) + "\t" + str(mds[0]) + "\t" + 
			str(mds[1]) + "\t" + str(mds[1] - mds[0] + 1) + "\t" + str(mds[2]) + "\n")
		
		# Close mds file
		mdsFile.close()
		
	# Update telomeres file
	telFile = open(Output_dir + '/Database_Input/tel.tsv', 'a')
	
	# Get number of telomeres
	tel_num = 0
	if left_Tel:
		tel_num += 1
	if right_Tel:
		tel_num += 1
	
	# Output to file
	telFile.write("\\N\t\\N\t" + contig + "\t" + str(mac_length) + "\t" + str(tel_num) + "\t")
	# Info about left telomere
	if left_Tel:
		telFile.write(str(left_Tel[0]+1) + "\t" + str(left_Tel[1]+1) + "\t" + str(left_Tel[1] - left_Tel[0] + 1) + "\t")
	else:
		telFile.write("\\N\t\\N\t0\t")
	# infor about right telomere
	if right_Tel:
		telFile.write(str(right_Tel[0]+1) + "\t" + str(right_Tel[1]+1) + "\t" + str(right_Tel[1] - right_Tel[0] + 1) + "\n")
	else:
		telFile.write("\\N\t\\N\t0\n")	
		
	telFile.close()

	# Update arrange table file
	arrFile = open(Output_dir + '/Database_Input/arr.tsv', 'a')
	
	# For each mic to mac map, output database entry
	for mic in MIC_to_HSP:
		# Get hsp list and declare variables
		hsp_list = sorted(MIC_to_HSP[mic], key=lambda x: int(x[7]) if int(x[7]) < int(x[8]) else int(x[8]))
		Arrangement = []
		nuc_shared = 0
		mismatch = 0
		dist_mds = set()
				
		# Iterate through hsp_list and get all needed information
		for hsp in hsp_list:
			Arrangement.append(-hsp[-1] if int(hsp[7]) > int(hsp[8]) else hsp[-1])
			nuc_shared += (int(hsp[3]) - int(hsp[4]))
			mismatch += int(hsp[4])
			dist_mds.add(hsp[-1])
		
		# Check if arrangement is scrambled
		is_scrambled = is_Scrambled(hsp_list, len(MDS_List), len(MDS_List)==len(dist_mds))
		# Put Arrangement into the canonical form
		Arrangement = toCanonicalForm(Arrangement, len(MDS_List))
		# Build arrangement string
		arrangement = ""
		for m in Arrangement[:-1]:
			if m > 0:
				arrangement += str(m) + ":0|" 
			else:
				arrangement += str(-m) + ":1|"
		if Arrangement[-1] > 0:
			arrangement += str(Arrangement[-1]) + ":0" 
		else:
			arrangement += str(-Arrangement[-1]) + ":1" 
			
		# Output arrangement table entry
		arrFile.write("\\N\t" + contig + "\t\\N\t" + str(mac_length) + "\t" + hsp_list[0][11] + "\t" + mic + "\t" + "\t\\N\t\\N\t" + str(nuc_shared) + "\t\\N\t" + 
					str(len(dist_mds)) + "\t" + str(mismatch) + "\t" + ("1" if is_scrambled else "0") + "\t" + arrangement + "\n")
	
	arrFile.close()
	
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function sorts hsp list

def sortHSP_List(MIC_maps):
		# Define sort function that will sort by:
		# 1) Higher Coverage, 2) Higher Length, 3) Higher Persent  identity match, 4) Lower Bitscore, 5) MIC, 6) the hsp start position in the MAC
		sort_func = key=lambda x: (float(x[11]), int(x[3]), float(x[2]), -float(x[10]), x[1], -int(x[5]))
		
		# Run sort
		MIC_maps.sort(key = sort_func, reverse=True)
	
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function runs telomeric regular expression, identifies telomere (if any) and stores telomeric sequences in the list
# Note: side = 5 is for 5' telomere and side = 3 is for 3' telomere

def identifyTelomere(reg_exp, seq, tel_seq, side):
	# Telomere to return
	tel_toReturn = None
	
	# Get the list of telomeric sequences
	telomeres = reg_exp.finditer(seq)
	tel_positions = sorted([(m.span()[0], m.span()[1]) for m in telomeres])
	# If this is a 5' telomeres
	if side == 5:
		# Go through each telomeric seq and build a telomere
		ind = 0
		for coord in tel_positions:
			if coord[1] - coord[0] >= Options['TelomereLength'] and coord[0] <= Options['TelomereEndLimit']:
				tel_toReturn = [coord[0], coord[1]]
				ind += 1
				break
			if coord[0] > Options['TelomereEndLimit']:
				break	
			ind += 1
			
		# check if left (5') telomeres can be extended
		if tel_toReturn and tel_positions:
			for tel in tel_positions[ind:]:
				# If two telomeric sequences are within tolerance error, merge them
				if tel[0] - tel_toReturn[1] > Options["TelomericErrorTolerance"]:
					break
				tel_toReturn[1] = tel[1]
	
	# If this is a 3' telomeres			
	elif side == 3:
		# Go through each telomeric seq and build a telomere
		ind = 0
		for coord in reversed(tel_positions):
			if coord[1] - coord[0] >= Options['TelomereLength'] and len(seq) - coord[1] <= Options['TelomereEndLimit']:
				tel_toReturn = [coord[0], coord[1]]
				ind += 1
				break
			if len(seq) - coord[1] > Options['TelomereEndLimit']:
				break	
			ind += 1
		# check if right (3') telomeres can be extended
		if tel_toReturn and tel_positions:
			for i in range(len(tel_positions) - ind - 1, -1, -1):
				tel = tel_positions[i]
				# If two telomeric sequences are within tolerance error, merge them
				if tel_toReturn[0] - tel[1] > Options["TelomericErrorTolerance"]:
					break
				tel_toReturn[0] = tel[0]
	
	# We have an error
	else:
		print("Error in identifyTelomere function, side = ", side, " while allowed values are 3 and 5")
		sys.exit()
	
	# Remove too short, non-telomeric sequences
	if tel_toReturn:
		tel_positions = [x for x in tel_positions if (x[1] - x[0] + 1) >= Options['TelomereLength'] or (tel_toReturn[0] <= x[0] and tel_toReturn[1] >= x[1])]
	else:
		tel_positions = [x for x in tel_positions if (x[1] - x[0] + 1) >= Options['TelomereLength']]
	# Append telomeric positions to tel_seq list
	tel_seq += tel_positions
	
	return tel_toReturn
	
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function updates mac_mds.gff3 file

def updateGFF(contig, MDS_List, Output_dir):
	# Open gff file
	gff = open(Output_dir + "/GFF/mac_mds.gff3", "a")
		
	# Output every mds
	for mds in MDS_List:
		gff.write(contig + "\tMI-ASS\tmds\t" + str(mds[0]) + "\t" + str(mds[1]) + "\t" + ".\t.\t.\tID=mds{0:06d};".format(updateGFF.mdsID) + "Name=mds_" + str(mds[-1]) + ";Target=" + contig + "\n")
		updateGFF.mdsID += 1
	
	gff.close()
	
updateGFF.mdsID = Options['MDS_id_start']
	
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function identifies best MIC contigs from which MAC was mapped and identifies scrambling

def identify_MIC_patterns(MIC_maps, MDS_List, MIC_to_HSP, Output_dir):
	if not MIC_maps:
		return
	# Variables for updating statistics on scrambled and complete
	stat_Scrambled = False
	stat_Complete = False
	stat_CompletScrambled = False
	
	# map MIC contigs to number of distinct mdss it has
	cont_to_mds = {}
	
	for mic in MIC_to_HSP:
		mdsNUM = len(set([x[-1] for x in MIC_to_HSP[mic]]))
		cont_to_mds[mic] = mdsNUM
	
	MICs = list(MIC_to_HSP.keys())
	# Sort by:
	# 1) The biggest number of distinct MDSs MIC has and 2) The highest MIC coverage
	MICs.sort(key=lambda x: (cont_to_mds[x], float(MIC_to_HSP[x][0][11])), reverse=True)
	
	out = open(Output_dir + "/Scrambling/all.tsv", "a")
	for mic in MICs:
		next = sorted(MIC_to_HSP[mic], key=lambda x: int(x[7]) if int(x[7]) < int(x[8]) else int(x[8]))
		out.write(MIC_maps[0][0] + "\t" + mic + "\t" + "{")
		for hsp in next[:-1]:
			out.write(("-" if int(hsp[7]) > int(hsp[8]) else "") + str(hsp[-1]) + ",")
		out.write(("-" if int(next[-1][7]) > int(next[-1][8]) else "") + str(next[-1][-1]) + "}\t")		
		
		# Check if this is a complete mapping
		if cont_to_mds[mic] == len(MDS_List):
			stat_Complete = True
			out.write("Complete\t")
		else:
			out.write("Incomplete\t")
		
		# check if it is a scrambled contig
		if(is_Scrambled(next, len(MDS_List), cont_to_mds[mic] == len(MDS_List))):
			stat_Scrambled = True
			out.write("Scrambled\n")
			
			# Output scrambled MIC pattern
			scramb_out = open(Output_dir + "/Scrambling/scrambled.tsv", "a")
			scramb_out.write(MIC_maps[0][0] + "\t" + mic + "\t" + "{")
			for hsp in next[:-1]:
				scramb_out.write(("-" if int(hsp[7]) > int(hsp[8]) else "") + str(hsp[-1]) + ",")
			scramb_out.write(("-" if int(next[-1][7]) > int(next[-1][8]) else "") + str(next[-1][-1]) + "}\t")
			
			# Check if this is a complete pattern
			if cont_to_mds[mic] == len(MDS_List):
				stat_CompletScrambled = True
				scramb_out.write("Complete\n")
			else:
				scramb_out.write("Incomplete\n")
			scramb_out.close()
		else:
			out.write("Non-Scrambled\n")
	
	# Select the best MIC to MAC maps taken from the sorting procedure 
	best_mic = MICs[0]
	hsp_list = sorted(MIC_to_HSP[best_mic], key=lambda x: int(x[7]) if int(x[7]) < int(x[8]) else int(x[8]))
	# Process this contig and its hsp 
	process_MIC_MAC_map(hsp_list, cont_to_mds[best_mic] == len(MDS_List), len(MDS_List), Output_dir)
	
	# If this is a complete map, then check if there are any other good mic maps and process them
	if cont_to_mds[best_mic] == len(MDS_List):
		for mic in MICs[1:]:
			if cont_to_mds[mic] != len(MDS_List):
				break
			next = sorted(MIC_to_HSP[mic], key=lambda x: int(x[7]) if int(x[7]) < int(x[8]) else int(x[8]))
			process_MIC_MAC_map(next, cont_to_mds[mic] == len(MDS_List), len(MDS_List), Output_dir)
	# Else, check for other MICs that have similar MDS number and MAC coverage
	else:
		rest_mic = [x for x in MICs if cont_to_mds[x] == cont_to_mds[best_mic] and float(MIC_to_HSP[x][0][11]) == float(MIC_to_HSP[best_mic][0][11])]
		rest_mic.remove(best_mic)
		for mic in rest_mic:
			next = sorted(MIC_to_HSP[mic], key=lambda x: int(x[7]) if int(x[7]) < int(x[8]) else int(x[8]))
			process_MIC_MAC_map(next, cont_to_mds[mic] == len(MDS_List), len(MDS_List), Output_dir)
	
	# Update stats
	if stat_Complete:
		identify_MIC_patterns.scrambled += 1
	if stat_Scrambled:
		identify_MIC_patterns.complete += 1
	if stat_CompletScrambled:
		identify_MIC_patterns.complete_scrambled += 1
			
	out.close()

identify_MIC_patterns.scrambled = 0
identify_MIC_patterns.complete = 0
identify_MIC_patterns.complete_scrambled = 0
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function checks whether a given arrangement given by hsp list is scrambled
# Note: is_complete can be set to False when it is not known if a map is complete

def is_Scrambled(MIC, mdsNum, is_complete):
	# Get string of MIC pattern
	s = ""
	if is_complete:
		# Build string directly
		for hsp in MIC[:-1]:
			s += ("-" if int(hsp[7]) > int(hsp[8]) else "") + str(hsp[-1]) + ","
		s += ("-" if int(MIC[-1][7]) > int(MIC[-1][8]) else "") + str(MIC[-1][-1])
	else:
		# Get set of MDSs and turn it into sorted list
		present_mdss = sorted(list({int(x[-1]) for x in MIC}))
		# Get MDS index map
		ind = 1
		MDS_map = dict()
		for mds in present_mdss:
			MDS_map[mds] = ind
			ind += 1
		# Update mdsNum
		mdsNum = ind - 1
		
		# Build string
		for hsp in MIC[:-1]:
			s += ("-" if int(hsp[7]) > int(hsp[8]) else "") + str(MDS_map[hsp[-1]]) + ","
		s += ("-" if int(MIC[-1][7]) > int(MIC[-1][8]) else "") + str(MDS_map[MIC[-1][-1]])
		
	# Get regular expressions for non-scrambled patterns
	r1 = ""
	for i in range(1, mdsNum):
		r1 += str(i) + ",(-?[0-9]*,)*"
	r1 += str(mdsNum)
	
	r2 = ""
	for i in range(mdsNum, 1, -1):
		r2 += "-" + str(i) + ",(-?[0-9]*,)*"
	r2 += "-1"
	#print("Reg exp 1: ", r1)
	#print("Reg exp 2: ", r2)
	
	# Check for non-scrambled pattern 1
	r1_comp = re.compile(r1)
	if r1_comp.search(s):
		return False
		
	# Check for non-scrambled pattern 2
	r2_comp = re.compile(r2)
	if r2_comp.search(s):
		return False
	
	# If program have not returned, then it is a scrambled pattern
	return True
	
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function saves current hsp as best MIC to MAC map, brings MIC arrangement into the standard form and reduces it

def process_MIC_MAC_map(hsp_list, is_complete, mdsNum, Output_dir):
	# First, get MIC arrangement
	Arrangement_0 = []
	if is_complete:
		# Build arrangement directly
		for hsp in hsp_list:
			Arrangement_0.append(-hsp[-1] if int(hsp[7]) > int(hsp[8]) else hsp[-1])
	else:
		# Get set of MDSs and turn it into sorted list
		present_mdss = sorted(list({int(x[-1]) for x in hsp_list}))
		# Get MDS index map
		ind = 1
		MDS_map = dict()
		for mds in present_mdss:
			MDS_map[mds] = ind
			ind += 1
		# Update mdsNum
		mdsNum = ind - 1
		
		# Build arrangement
		for hsp in hsp_list:
			Arrangement_0.append(-MDS_map[hsp[-1]] if int(hsp[7]) > int(hsp[8]) else MDS_map[hsp[-1]])
	# Remove consecutive repeating letters (ex: 1, 2, 3, 3, 4, 5, 5, 5 - > 1, 2, 3, 4, 5)
	Arrangement = [Arrangement_0[0]]
	if len(Arrangement_0) > 1:
		prev = Arrangement_0[0]
		for m in Arrangement_0[1:]:
			if m != prev:
				Arrangement.append(m)
				prev = m
	
	# Get arrangement in the canonical form
	Arrangement = toCanonicalForm(Arrangement, mdsNum)
	
	# Get reduce arrangement
	reduced = [Arrangement[0]]
	if len(Arrangement) > 1:
		prev = Arrangement[0]
		for m in Arrangement[1:]:
			# If both positive and increasing, continue
			if m > 0 and prev > 0 and m == prev + 1:
				prev = m
				continue
			# If both negative and decreasing, continue
			elif m < 0 and prev < 0 and m == prev + 1:
				prev = m
				continue
			# Else, append it to the reduced list
			reduced.append(m)
			prev = m
	
	# Get index map
	ind = 1
	Ind_map = dict()
	for mds in sorted(list({abs(x) for x in reduced})):
		Ind_map[mds] = ind
		ind += 1
	# Update mdsNum
	mdsNum = ind - 1
	
	# Map mdss to for reduced arrangement
	Reduced = []
	for mds in reduced:
		Reduced.append(-Ind_map[abs(mds)] if mds < 0 else Ind_map[abs(mds)])
	
	# Put reduced arrangement into the canonical form
	Reduced = toCanonicalForm(Reduced, mdsNum)
	
	# Output result
	out = open(Output_dir + "/Scrambling/maps.tsv", "a")
	out.write(hsp_list[0][0] + "\t" + hsp_list[0][1] + "\t" + "{" + arrangementToString(Arrangement) + "}\t{" + arrangementToString(Reduced) + "}\n")
	out.close()
	
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function counts number of inverted MDSs in the arrangement

def getNumber_Inv_MDS(Arrangement):
	count = 0
	for m in Arrangement:
		if m < 0:
			count += 1
	return count

#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function puts arrangement to string

def	arrangementToString(Arrangement):
	s = ""
	for m in Arrangement[:-1]:
		s += str(m) + ", "
	s += str(Arrangement[-1])
	return s

#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function returns the position of the first inversion in the arrangement

def firstInv(Arrangement):
	for i in range(0, len(Arrangement)):
		if Arrangement[i] < 0:
			return i+1
	# No inversions
	return 0

#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function takes MIC arrangement return the canonical form of it

def toCanonicalForm(Arrangement, mdsNum):
	# Put arrangement in the right order:
	# 1) minimal number of inversions, 2) lowest lexicographical order, 3) delay first inversion
	# Get the other 3 arrangements first
	Arrangement_I = []
	for i in range(len(Arrangement) - 1, -1, -1):
		m = Arrangement[i]
		Arrangement_I.append(-1 * m)
		
	Arrangement_A = []
	for m in Arrangement:
		if m > 0:
			Arrangement_A.append(-1 * (mdsNum + 1 - m))
		else:
			Arrangement_A.append(mdsNum + 1 - abs(m))
			
	Arrangement_AI = []
	for i in range(len(Arrangement_A) - 1, -1, -1):
		m = Arrangement_A[i]
		Arrangement_AI.append(-1 * m)
	# Put all arrangements in the list and sort it by above criterias
	Arrangement_List = [Arrangement, Arrangement_I, Arrangement_A, Arrangement_AI]
	Arrangement_List.sort(key=lambda x: (getNumber_Inv_MDS(x), arrangementToString(x), -firstInv(x)))
	return Arrangement_List[0]

#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function takes MIC maps and MDS list and removes noise hsps from MIC_maps

def removeNoise(MIC_maps, MDS_List):
	if not MDS_List or not MIC_maps:
		return
	# First make a dictionary of MIC to its hsps
	MIC_to_HSP = {}
	for hsp in MIC_maps:
		if hsp[1] in MIC_to_HSP:
			MIC_to_HSP[hsp[1]].append(hsp)
		else:
			MIC_to_HSP[hsp[1]] = [hsp]
			
	# List of MICs to remove
	MICtoRemove = []
	
	# Mark MICS that does not meet MDS threshold requirement
	# 1) Check if there are any hsp that satisfies min hsp to mds ratio
	# 2) Check if mds number percentage is satisfied
	mdsNum = math.log(len([x for x in MDS_List if x[2] == 0]))
	for mic in MIC_to_HSP:
		hsp_list = MIC_to_HSP[mic]
		is_Good = False
		# Check for condition 1)
		for hsp in hsp_list:
			mds = MDS_List[hsp[-1] - 1]
			if (mds[1] -  mds[0] + 1)/(float(hsp[3])) >= Options['Min_hsp_to_mds_ratio']:
				is_Good = True
				break
		if is_Good:
			continue
		
		# Check for condition 2)
		distMDS = len({x[-1] for x in hsp_list})	
		if math.log(distMDS)/mdsNum >= Options['Min_mds_num_percentage']:
			continue
			
		# Else, add MIC as the one for removal
		MICtoRemove.append(mic)
		
	# Filter MIC contigs that are bad
	MIC_maps[:] = [x for x in MIC_maps if x[1] not in MICtoRemove]
			
			
			
			
			
			
			
			
			
			
			
			
			





