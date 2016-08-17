import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import os

#build a connectivity matrix for each run
def buildMatrix(matrixFile):
	inputH=[]
	with open(matrixFile) as input:
		for line in input:
			line=line.rstrip('\n')
			line=line.split('\t')
			converted=[float(x) for x in line]
			inputH.append(converted)
	
	inputH=np.array(inputH)
	sampleNum=len(inputH[0])
	connectivityMat=np.zeros((sampleNum, sampleNum), dtype=np.int)

	#fills in connectivity matrix
	for connectivityX in range(0, sampleNum):
		for connectivityY in range(0, sampleNum):
			#return cluster location of the maximum metagene factor for sample x vs sample y
			locationOfMax1=list(inputH[:,connectivityX]).index(max(list(inputH[:,connectivityX])))
			locationOfMax2=list(inputH[:,connectivityY]).index(max(list(inputH[:,connectivityY])))
			#if the maximum metagene factor for sample x and sample y are in the same cluster location, they are connected
			if locationOfMax1==locationOfMax2:
				connectivityMat[connectivityX, connectivityY]=1

	#outputs connectivity matrix in tab-delimited format
	outputMatrix=open(str(matrixPath)+'final_connectivity_matrix'+str(args.matrixFile.split("/")[-1]), 'w')
	for i in range(0, sampleNum):
		for n in range(0, sampleNum-1):
			outputMatrix.write(str(connectivityMat[i][n])+'\t')
		outputMatrix.write(str(connectivityMat[i][sampleNum-1])+'\n')

	return connectivityMat

def visualize_connectivity(connectivityMat, sampleName):
	if sampleName=='noXLabels':
		fig1, ax1=plt.subplots()
		heatmap=ax1.pcolor(connectivityMat, cmap=plt.cm.jet)
		cbar = plt.colorbar(heatmap)
		plt.ylabel('sample ID')
		plt.xlabel('sample ID')
		plt.savefig(str(visPath)+'final_connectivity_matrix'+str(args.matrixFile.split("/")[-1])+'.png')
	#if column names were provided by user, heatmap is labeled
	else:
		colNames=[]
		with open(sampleName) as input:
			for line in input:
				colNames.append(line.rstrip('\n'))
		fig1, ax1=plt.subplots()
		heatmap=ax1.pcolor(connectivityMat, cmap=plt.cm.jet)
		cbar = plt.colorbar(heatmap)
		ax1.set_yticks(np.arange(connectivityMat.shape[1])+0.5, minor=False)
		ax1.set_yticklabels(colNames, minor=False)
		ax1.set_xticks(np.arange(connectivityMat.shape[1])+0.5, minor=False)
		ax1.set_xticklabels(colNames, minor=False)
		plt.ylabel('sample ID')
		plt.xlabel('sample ID')
		plt.savefig(str(visPath)+'final_connectivity_matrix'+str(args.matrixFile.split("/")[-1])+'.png')

if __name__=='__main__':
	parser=argparse.ArgumentParser("parses information to build connectivity matrix")
	parser.add_argument('-input', required=True, dest='matrixFile', help='Full path to tab-delimited "H matrix" file')
	parser.add_argument('--colNames', default='noXLabels', dest='colNames', type=str, help='full path to file of sample names in order of matrix, one name per line')
	parser.add_argument('--output', default=os.getcwd(), dest='outPath', type=str, help='full path to output directory')
	args=parser.parse_args()
	
	# path to output directories
	visPath=str(args.outPath)+'connectivity_visualization/'
	matrixPath=str(args.outPath)+'connectivity_matrix/'

	#check if output directory is already made, if not creates it
	if os.path.isdir(visPath) == False:
		os.mkdir(visPath)
	if os.path.isdir(matrixPath) == False:
		os.mkdir(matrixPath)

	#input is a file of predicted matrix H (kxm, k=clusters/metagene expression, m=samples)
	connectivityMat=buildMatrix(matrixFile=args.matrixFile);
	visualize_connectivity(connectivityMat=connectivityMat, sampleName=args.colNames);