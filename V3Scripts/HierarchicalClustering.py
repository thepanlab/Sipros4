
#!/usr/bin/python


import sys, getopt, warnings, os, re




def mergeClosestClusters (CandidateCluster, distanceThreshold) :
	bReVal = True
	firstBestCulsterId = -1
	secondBestCulsterId = -1
	#print CandidateCluster
	for i in range(len(CandidateCluster) - 1) :
		averageFirst  = sum(CandidateCluster[i])/float(len(CandidateCluster[i]))
		averageSecond = sum(CandidateCluster[i+1])/float(len(CandidateCluster[i+1]))
		if (averageFirst > averageSecond) :
			print "wrong rank!"
			sys.exit(0) 
		currentDistance = averageSecond - averageFirst
		if ( currentDistance  <= distanceThreshold) :
			if ((firstBestCulsterId == -1) or (secondBestCulsterId == -1)): # first pair of good clusters
				minDistance = currentDistance
				firstBestCulsterId = i
				secondBestCulsterId = i+1
			elif ( currentDistance  < minDistance)   : # two bettter clusters
				minDistance = currentDistance
                                firstBestCulsterId = i
                                secondBestCulsterId = i+1
		#	print minDistance, currentDistance
	if ((firstBestCulsterId != -1) and (secondBestCulsterId != -1)) :
		#merge two clusters
		mergedCluster = CandidateCluster [firstBestCulsterId] + CandidateCluster [secondBestCulsterId]
		del CandidateCluster[firstBestCulsterId]
		del CandidateCluster[firstBestCulsterId]
		CandidateCluster.insert(firstBestCulsterId, mergedCluster)
	else :
		bReVal = False
		
	return bReVal


def hierarchicalClustering (ldCandidatePct, distanceThreshold) :
	ldCandidatePct.sort()
	CandidateCluster = []
	if (len(ldCandidatePct) == 1) :
		CandidateCluster.append([ldCandidatePct[0]])
	elif (len(ldCandidatePct) > 1) :
		# each cluster has one candidate 
		for i in range(len(ldCandidatePct)) : 
			CandidateCluster.append([ldCandidatePct[i]])
		#cluters merge 
		bMerge =  mergeClosestClusters(CandidateCluster, distanceThreshold)
		while (bMerge) :
		#	print CandidateCluster
			bMerge =  mergeClosestClusters(CandidateCluster, distanceThreshold)
	return CandidateCluster
