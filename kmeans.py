import argparse
import numpy as np

def initializeCentroids(inputMatrix, k):
	inputH=[]
	with open(args.matrixFile) as input:
		for line in input:
			line=line.rstrip('\n')
			line=line.split('\t')
			converted=[float(line[x]) for x in range(0, len(line))]
			inputH.append(converted)
	
	inputH=np.array[inputH]
	#randomly chooses k samples for centroid initiation point
	startingCentroids=random.sample(xrange(1, len(inputH[0])), k)
	#create centroid and cluster groups, key=cluster number, value[0]=values of centroid (feature values), value[1]=list of patients in cluster
	clusters={'clusters_'+str(i):[inputH[,startingCentroids[i]], []] for i in range(0, k)}
	
	return inputH, clusters


def calcEuclideanDist(inputH, clusters):
	for i in range(0, len(inputH[0])):
		patientFeatures=inputH[,i]
		#stores the Euclidean distance for each cluster for an individual patient
		eucDistCalcs=[]
		for j in clusters:
			temp=0
			#index 0 always contains centroid values
			for centroidValues in range(0, len(j[0])):
				temp=temp+math.pow(j[0][centroidValues]-patientFeatures[centroidValues], 2)
			eucDistCalcs.append((j, math.sqrt(temp)))
		#returns tuple with smallest Euc distance
		clusterAssignment=min(eucDistCalcs, key=lamda x:x[1])
		#TO DO:  if more than one cluster assignment is returned, need to figure out which cluster to assign to
		#add patientID to cluster
		clusters[clusterAssignment[0]]=clusters[clusterAssignment[0]]+[str(i)]

	return clusters

def updateCentroids(clusters):
	pass;

if __name__=='__main__':
	parser=argparse.ArgumentParser(description='kmeans clustering of NMF matrix H')
	parser.add_argument('-input', required=True, dest='matrixFile', help='Full path to tab-delimited "H matrix" file')
	parser.add_argument('-kclusters', default='2', dest='kclusters', type=int, help='[INT] Number of subtypes or clusters to expect, must be smaller than m columns and n rows of input data')
	parser.add_argument('-iterations', default='100', dest='iterations', type=int, help='[INT] Number of iterations requried for convergence')
	args=args=parser.parse_args()
	inputH, clusters=initializeCentroids(inputMatrix=args.matrixFile, k=args.kclusters);
	
	for i in range(0, args.kclusters):	
		clusters=calcEuclideanDist(inputH, clusters);