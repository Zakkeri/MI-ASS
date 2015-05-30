# Python file with all essential functions used during CiliateAnnotation program execution

#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function goes through the list of high scoring pairs that are associated with the MIC and constructs MDSs for the MAC

def getMDS_List(MIC_maps):
	MDS_List = list()
	
	#print('Before sorting: ', MIC_maps)
	# Sort MIC maps by coverage
	MIC_maps.sort(key=lambda x: float(x[0][11]), reverse=True)
	#print('After sorting: ',MIC_maps)
	#print('Working with ', str(MIC_maps[0][0]), '\n')
	# Build list of MDSs
	for mic in MIC_maps:
		for hsp in mic:
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


#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function checks whether there are any gaps in the MAC annotation

def addGaps(MDS_List, MAC_start, MAC_end):

	# Get the list of covered MAC Interval(s)
	MAC_Interval = list()
	for mds in sorted(MDS_List, key = lambda x: x[0]):
		if not MAC_Interval:
			MAC_Interval.append([mds[0], mds[1], mds[2]])
		else:
			if MAC_Interval[-1][1] >= mds[0]:
				MAC_Interval[-1][1] = max(mds[1], MAC_Interval[-1][1])
			else:
				MAC_Interval.append([mds[0], mds[1], mds[2]])
	
	# If we have more than one interval, then there are gaps we need to add to the annotation
	if len(MAC_Interval) > 1:
		prev = MAC_Interval[0]
		for interv in MAC_Interval[1:]:
			MDS_List.append([prev[1], interv[0],0])
			prev = interv
			
	# Check for gaps at the begining of MAC and at the end of MAC
	if MAC_Interval[0][0] - MAC_start > 0:
		MDS_List.append([MAC_start, MAC_Interval[0][0]], 0)
	if MAC_end - MAC_Interval[-1][1] > 0:
		MDS_List.append([MAC_Interval[0][1]], MAC_end, 0)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
