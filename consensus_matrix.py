import argparse
import numpy as np
import pandas as pd
import seaborn as sns

def readMatrices(inputFile):
	#stores connectivity matrix for each run key=runID:value=connectivity matrix
	connectivityMatrices={}
	with open(inputFile) as name:
		for filename in name:
			oneMatrix=[]
			with open(filename.rstrip('\n')) as input:
				for line in input:
					line=line.rstrip('\n')
					line=line.split('\t')
					converted=[float(line[x]) for x in range(0, len(line))]
					oneMatrix.append(converted)
				#stores the connectivity matrix for each run	
				connectivityMatrices[filename]=np.array(oneMatrix)
	
	dimOfConsensus=len(oneMatrix[0])
	#creates consensus matrix of all zeroes
	consensusMat=np.zeros((dimOfConsensus, dimOfConsensus), dtype=np.float)
	
	return connectivityMatrices, consensusMat		

def buildConsensus(connectivityMatrices, consensusMat):

	#iterates through all matrix indices so takes the average of all numbers of connectivity matrix in that
	#position to create consensus matrix
	for i in range(0, len(consensusMat[0])):
		for j in range(0, len(consensusMat[0])):
			total=sum([connectivityMatrices[key][i][j] for key in connectivityMatrices])
			average=total/float(len(connectivityMatrices))
			consensusMat[i][j]=average

	return consensusMat

def visualizeConsensus(consensusMat, colNames):
	if colNames=='noXLabels':
		#put concensus matrix into dataframe to build hierarchical clustermap		
		dataframe=pd.DataFrame(data=consensusMat)
		#clusters by columns and rows and annotates probablility a particular sample clusters together
		#cluster distance is meausred by average Euclidean Distance in seaborn for hierarchical clustering
		concensusClustered=sns.clustermap(dataframe, col_cluster=True, row_cluster=True, annot=True)
		sns.plt.show()
	
	else:
		#assigns sample names to consensus matrix
		sampleNames=[]
		with open(colNames) as input:
			for line in input:
				sampleNames.append(line.rstrip('\n'))
		#put concensus matrix into dataframe to build hierarchical clustermap		
		dataframe=pd.DataFrame(data=consensusMat, index=sampleNames, columns=sampleNames)
		#clusters by columns and rows and annotates probablility a particular sample clusters together
		#cluster distance is meausred by average Euclidean Distance in seaborn for hierarchical clustering
		concensusClustered=sns.clustermap(dataframe, col_cluster=True, row_cluster=True, annot=True)
		sns.plt.show()

if __name__=='__main__':
	parser=argparse.ArgumentParser(description='builds consensus matrix given set of connectivity matrices')
	parser.add_argument('-input', required=True, dest='listOfMatrices', help='List of connectivity matrices')
	parser.add_argument('-colNames', default='noXLabels', dest='colNames', type=str)
	args=parser.parse_args()
	
	connectivityMatrices, consensusMat=readMatrices(inputFile=args.listOfMatrices);
	consensusMat=buildConsensus(connectivityMatrices, consensusMat)
	visualizeConsensus(consensusMat, colNames=args.colNames)
