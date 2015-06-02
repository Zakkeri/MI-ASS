# Python file with all essential functions used during CiliateAnnotation program execution
import subprocess
import sys
from settings import *

#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function runs rough BLAST and returns the hsp result list

def run_Rough_BLAST(Output_dir, contig):
	# Set up parameters
	dust = "yes" if Options['RoughBlastDust'] else "no"
	ungapped = " -ungapped " if Options['RoughBlastUngapped'] else ""
	maskLowercase = " -lcase_masking " if Options['BlastMaskLowercase'] else ""
	
	# Run BLAST command
	rough_out = subprocess.check_output("blastn -task " + Options['RoughBlastTask'] + " -word_size " + str(Options['RoughBlastWordSize']) + " -max_hsps 0 " +
	"-max_target_seqs 10000 -dust " + dust + ungapped + maskLowercase + "-query " + Output_dir + "/hsp/rough/masked_" + str(contig) + ".fa -db " +
	Output_dir + "/blast/mic -num_threads " + str(Options['ThreadCount']) + 
	" -outfmt \"10 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qcovs\"")
	
	# Filter empty rows
	roughVal = [x.rstrip() for x in rough_out.decode(sys.stdout.encoding).split('\n') if x != ""]
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
	dust = "yes" if Options['FineBlastDust'] else "no"
	ungapped = " -ungapped " if Options['FineBlastUngapped'] else ""
	maskLowercase = " -lcase_masking " if Options['BlastMaskLowercase'] else ""
	
	fine_out = subprocess.check_output("blastn -task " + Options['FineBlastTask'] + " -word_size " + str(Options['FineBlastWordSize']) + " -max_hsps 0 " + 
	"-max_target_seqs 10000 -dust " + dust + ungapped + maskLowercase + "-query " + Output_dir + "/hsp/rough/masked_" + str(contig) + ".fa " +
	"-db " + Output_dir + "/blast/mic " +
	"-outfmt \"10 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qcovs\"" + " -num_threads " + str(Options['ThreadCount']))
	
	# Filter empty rows
	fineVal = [x.rstrip() for x in fine_out.decode(sys.stdout.encoding).split('\n') if x != ""]
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

# This function goes through the list of high scoring pairs that are associated with the MIC and constructs MDSs for the MAC

def get_Rough_MDS_List(MIC_maps):
	MDS_List = list()
	
	#print('Before sorting: ', MIC_maps)
	# Sort MIC maps by coverage
	#MIC_maps.sort(key=lambda x: float(x[0][11]), reverse=True)
	#print('After sorting: ',MIC_maps)
	#print('Working with ', str(MIC_maps[0][0]), '\n')
	# Build list of MDSs
	for hsp in MIC_maps:
		mds_toAdd = [int(hsp[5]), int(hsp[6]),1]
		#print("Considering: ", mds_toAdd)
		# Check if it is a subset of some MDS and skip it if it does
		if [x for x in MDS_List if mds_toAdd[0] >= x[0] and mds_toAdd[1] <= x[1]]:
			continue
			
		# Get list of hsp that overlap with hsp that we are trying to add
		overlap = [x for x in MDS_List if (x[0] > mds_toAdd[0] and x[0] < mds_toAdd[1]) or (x[1] > mds_toAdd[0] and x[1] < mds_toAdd[1])]
		toAdd = False
		
		# If overlap is empty, then add MDS
		if not overlap:
			#print("No overlap, just add")
			toAdd = True
		# Go through current MDSs and see if any can be made longer
		else:
			for x in sorted(overlap, key=lambda x: x[0]):
				# Check if two MDSs can be merged
				if (mds_toAdd[0] + mds_toAdd[1])/2 in range(x[0], x[1]):
					#print("Can merge two MDSs")
					mds_toAdd[0] = min(mds_toAdd[0], x[0])
					mds_toAdd[1] = max(mds_toAdd[1], x[1])
					#print("Removing: ", x)
					MDS_List.remove(x)
					toAdd = True
				# Check if non-overlapping hsp portion overlaps with some other MDS and if it doesn't, then add it
				# Overlap is on the left
				elif mds_toAdd[0] < x[0]:
					#print("Overlap on the left")
					sub_overlap = [y for y in MDS_List if (y[1] > mds_toAdd[0]) and (y[1] < x[1])]
					if not sub_overlap:
						toAdd = True
				# Overlap is on the right
				else:
					#print("Overlap on the right")
					sub_overlap = [y for y in MDS_List if (y[0] > x[0] and y[0] < mds_toAdd[1])]
					if not sub_overlap:
						toAdd = True
						
		if toAdd:
			#print("Adding: ", mds_toAdd)
			MDS_List.append(mds_toAdd)
	
	return MDS_List
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function takes fine BLAST result and current MDS list and tries to fill the gaps and improve annotation

def improveAnnotation(Fine_BLAST, MDS_List, MAC_Coverage, MAC_start, MAC_end):
	# Improve annotation iteratively by trying to fill the gaps untill no gaps, or no changes
	is_Change = True
	while is_Change:
		is_Change = False
		
		# Build a list of gaps
		gaps = list()
		if MAC_Coverage[0][0] - MAC_start > 0:
			gaps.append([MAC_start, MAC_Coverage[0][0]])
		prev = MAC_Coverage[0]
		for interval in MAC_Coverage[1:]:
			gaps.append([prev[1], interval[0]])
			prev = interval
		if MAC_end - MAC_Coverage[-1][1] > 0:
			gaps.append([MAC_Coverage[-1][1], MAC_end])
		
		# Go through each gap and try to fill it with the hsp from fine BLAST output
		for gap in gaps:
			# Construct list of hsps that overlap with the gap
			gap_overlaps = [x for x in Fine_BLAST if gap[0] < int(x[6]) and gap[1] > int(x[5])]
			if not gap_overlaps:
				continue
		
			# Sort hsps by 1)covers most of the gap, 2) has higher coverage, 3) has lower bitscore
			reduce_func = lambda a: (min(int(a[6]), gap[1]) - max(int(a[5]), gap[0]),  float(a[11]), -float(a[10]))
			gap_overlaps.sort(key = reduce_func, reverse=True)
		
			# Go through the list and get the "best" fitting hsp
			hsp = gap_overlaps[0]
			mds_toAdd = [int(hsp[5]), int(hsp[6]),1] 
			# Add matched hsp to the list of MDSs
			MDS_List.append(mds_toAdd)
			is_Change = True
		
		# Sort the MDS List
		MDS_List.sort(key=lambda x: x[0])
		
		# Recalculate MAC_Coverage and check if we are done
		MAC_Coverage = getCovering_Intervals(MDS_List)
		if len(MAC_Coverage) == 1 and MAC_Coverage[0][0] - MAC_start <= 0 and MAC_end - MAC_Coverage[-1][1] <= 0:
			break
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This function takes a list of MDSs (sorted by the MDS begining coordinate) and returns the intervals of the MAC covering

def getCovering_Intervals(MDS_List):
	# Construct intervals by using MDS List and checking for gaps between consecutive MDSs
	MAC_Interval = list()
	for mds in sorted(MDS_List, key = lambda x: x[0]):
		if not MAC_Interval:
			MAC_Interval.append([mds[0], mds[1]])
		else:
			if MAC_Interval[-1][1] >= mds[0]:
				MAC_Interval[-1][1] = max(mds[1], MAC_Interval[-1][1])
			else:
				MAC_Interval.append([mds[0], mds[1]])
	
	return MAC_Interval
				
#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function checks whether there are any gaps in the MAC annotation

def addGaps(MDS_List, MAC_start, MAC_end):
	# If MDS List is empty, then return the whole MAC interval as a gap
	if not MDS_List:
		MDS_List.append([MAC_start, MAC_end, 0])
		return
	
	# Get the list of covered MAC Interval(s)
	MAC_Interval = getCovering_Intervals(MDS_List)
		
	# If we have more than one interval, then there are gaps we need to add to the annotation
	if len(MAC_Interval) > 1:
		prev = MAC_Interval[0]
		for interv in MAC_Interval[1:]:
			MDS_List.append([prev[1], interv[0],0])
			prev = interv
			
	# Check for gaps at the begining of MAC and at the end of MAC
	if MAC_Interval[0][0] - MAC_start > 0:
		MDS_List.append([MAC_start, MAC_Interval[0][0], 0])
	if MAC_end - MAC_Interval[-1][1] > 0:
		MDS_List.append([MAC_Interval[0][1], MAC_end, 0])
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
