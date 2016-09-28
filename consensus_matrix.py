import argparse
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def readMatrices(inputFile):
	#stores connectivity matrix for each run key=runID:value=connectivity matrix
	connectivityMatrices={}
	with open(inputFile) as name:
		#skips header line
		header=next(name)
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
			#extracts clusters number for file name creations
			clusters=filename.split('=')[1][0]

	dimOfConsensus=len(oneMatrix[0])
	#creates consensus matrix of all zeroes
	consensusMat=np.zeros((dimOfConsensus, dimOfConsensus), dtype=np.float)
	
	return connectivityMatrices, consensusMat, clusters		

def buildConsensus(connectivityMatrices, consensusMat):

	#iterates through all matrix indices so takes the average of all numbers of connectivity matrix in that
	#position to create consensus matrix
	for i in range(0, len(consensusMat[0])):
		for j in range(0, len(consensusMat[0])):
			total=sum([connectivityMatrices[key][i][j] for key in connectivityMatrices])
			average=total/float(len(connectivityMatrices))
			consensusMat[i][j]=average

	return consensusMat

def visualizeConsensus(consensusMat, connectivityMatrices, clusters, colNames):
	if colNames=='noXLabels':
		#put concensus matrix into dataframe to build hierarchical clustermap		
		dataframe=pd.DataFrame(data=consensusMat)
		#clusters by columns and rows and annotates probablility a particular sample clusters together
		#cluster distance is meausred by average Euclidean Distance in seaborn for hierarchical clustering
		consensusClustered=sns.clustermap(dataframe, col_cluster=True, row_cluster=True, annot=True)
		consensusClustered.savefig(str(matrixPath)+'consensus_Matrix_over_'+str(len(connectivityMatrices))+'_runs_at_k='+str(clusters)+'.png')
	
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
		consensusClustered=sns.clustermap(dataframe, col_cluster=True, row_cluster=True, annot=True)
		consensusClustered_non_annt=sns.clustermap(dataframe, col_cluster=True, row_cluster=True, annot=False)
		plt.setp(consensusClustered.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
		plt.setp(consensusClustered_non_annt.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
		plt.setp(consensusClustered.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
		plt.setp(consensusClustered_non_annt.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
		consensusClustered.savefig(str(matrixPath)+'consensus_Matrix_over_'+str(len(connectivityMatrices))+'_runs_at_k='+str(clusters)+'.png')
		consensusClustered_non_annt.savefig(str(matrixPath)+'non_annotated_consensus_Matrix_over_'+str(len(connectivityMatrices))+'_runs_at_k='+str(clusters)+'.png')

if __name__=='__main__':
	parser=argparse.ArgumentParser(description='builds consensus matrix given set of connectivity matrices')
	parser.add_argument('-input', required=True, dest='listOfMatrices', help='List of connectivity matrices')
	parser.add_argument('--colNames', default='noXLabels', dest='colNames', type=str)
	parser.add_argument('--output', default=os.getcwd(), dest='outPath', type=str, help='full path to output directory')
	args=parser.parse_args()
	
	# path to output directories
	matrixPath=str(args.outPath)+'consensus/'

	#check if output directory is already made, if not creates it
	if os.path.isdir(matrixPath) == False:
		os.mkdir(matrixPath)


	connectivityMatrices, consensusMat, clusters=readMatrices(inputFile=args.listOfMatrices);
	consensusMat=buildConsensus(connectivityMatrices, consensusMat)
	visualizeConsensus(consensusMat, connectivityMatrices, clusters, colNames=args.colNames)
