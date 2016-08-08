from time import gmtime, strftime
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os

def matrixInitialization(matrixFile, k):
	#reads in matrix data
	observed=[]
	with open(matrixFile) as input:
		for line in input:
			line=line.rstrip('\n')
			line=line.split('\t')
			converted=[float(line[x]) for x in range(0, len(line))]
			observed.append(converted)

	observed=np.array(observed)
	nrows=len(observed)
	mcols=len(observed[0])

	#randomly initialize matrix W and H to small non-negative numbers
	W=np.random.rand(nrows, k)
	H=np.random.rand(k, mcols)
	return observed, W, H

def updateRules(observed, W, H):
	'''-----------first update matrix W-------------
	NOTE: output of the new H in the input of W'''
	divide_W=np.divide(np.dot(observed,H.T), np.dot(np.dot(W,H),H.T))
	#new matrix W, elemental multiply, not dot product
	W=np.multiply(W, divide_W)

	'''----------then update matrix H--------------
	NOTE:  output of new W is input of H'''
	divide_H=np.divide(np.dot(W.T, observed), np.dot(np.dot(W.T, W), H))
	#new matrix H, elemental multiply, not dot product
	H=np.multiply(H, divide_H)
	#predicted value of input, **should be very close to input matrix**
	predicted=np.dot(W, H)
	return W, H, predicted

def costFunction(observed, predicted):
	#subtracts matricies observed from predicted, takes the sqaure of each element in subtracted matrix
	#and sums all elements together to get squared Euclidean distance
	diff=np.subtract(observed, predicted)
	squaredDiff=np.square(diff)
	squaredEucDist=np.sum(squaredDiff)
	return squaredEucDist

def visualizeConvergenceAccuracy(qualityApprox, iterConverge):
	#a line graph of the minimized costFunction across all iterations
	fig1, ax1=plt.subplots()
	plt.plot(iterConverge, qualityApprox)
	plt.title('squared Euclidean distance between observed vs predicted')
	plt.axis([0, max(iterConverge), 0, max(qualityApprox)+50])   
	plt.ylabel('squared Euclidean Distance')
	plt.xlabel('iterations')
	plt.savefig(str(args.outPath)+'sqEucDist_over_all_iterations_k='+str(args.kclusters)+'_'+str(uniqueName)+'.png')

	#zoomed in plot of the first 10 iterations
	plt.plot(iterConverge, qualityApprox)
	plt.title('squared Euclidean distance between observed vs predicted')
	plt.axis([0, 10, 0, max(qualityApprox)+50])
	plt.ylabel('squared Euclidean Distance')
	plt.xlabel('iterations')
	plt.savefig(str(args.outPath)+'sqEucDist_over_first_10_iterations_k='+str(args.kclusters)+'_'+str(uniqueName)+'.png')

def matrix_visualization(W, H, yAxisNames, xAxisNames):
	#construct heat map of matrix W
	def matrix_W():
		#if no names are provided for the rows it makes a heatmap with no labels
		if yAxisNames=='noYLabels':
			fig2, ax2=plt.subplots()
			heatmap=ax2.pcolor(W, cmap=plt.cm.RdYlGn)
			cbar = plt.colorbar(heatmap)
			plt.ylabel('genes')
			plt.xlabel('gene expression profiles (k clusters)')
			plt.savefig(str(args.outPath)+'matrixW_visualization_k='+str(args.kclusters)+'_'+str(uniqueName)+'.png')
		#if row names were provided by user, heatmap is labeled
		else:
			rowNames=[]
			with open(yAxisNames) as input:
				for line in input:
					rowNames.append(line.rstrip('\n'))
			fig2, ax2=plt.subplots()
			heatmap=ax2.pcolor(W, cmap=plt.cm.RdYlGn)
			cbar = plt.colorbar(heatmap)
			ax2.get_xaxis().set_visible(False)
			ax2.set_yticks(np.arange(W.shape[0])+0.5, minor=False)
			ax2.set_yticklabels(rowNames, minor=False)
			plt.ylabel('genes')
			plt.xlabel('gene expression profiles (k clusters)')
			plt.savefig(str(args.outPath)+'matrixW_visualization_k='+str(args.kclusters)+'_'+str(uniqueName)+'.png')
	
	#construct heatmap fo matrix H
	def matrix_H():
	#if no names are provided for the columns it makes a heatmap with no labels
		if xAxisNames=='noXLabels':
			fig3, ax3=plt.subplots()
			heatmap=ax3.pcolor(H, cmap=plt.cm.RdYlGn)
			cbar = plt.colorbar(heatmap)
			plt.ylabel('gene expression profiles (k clusters)')
			plt.xlabel('sample ID')
			plt.savefig(str(args.outPath)+'matrixH_visualization_k='+str(args.kclusters)+'_'+str(uniqueName)+'.png')
	#if column names were provided by user, heatmap is labeled
		else:
			colNames=[]
			with open(xAxisNames) as input:
				for line in input:
					colNames.append(line.rstrip('\n'))
			fig3, ax3=plt.subplots()
			heatmap=ax3.pcolor(H, cmap=plt.cm.RdYlGn)
			cbar = plt.colorbar(heatmap)
			ax3.get_yaxis().set_visible(False)
			ax3.set_xticks(np.arange(H.shape[1])+0.5, minor=False)
			ax3.set_xticklabels(colNames, minor=False)
			plt.ylabel('gene expression profiles (k clusters)')
			plt.xlabel('sample ID')
			plt.savefig(str(args.outPath)+'matrixH_visualization_k='+str(args.kclusters)+'_'+str(uniqueName)+'.png')

	matrix_W();
	matrix_H();


if __name__=='__main__':
	#arguments for classifier
	parser=argparse.ArgumentParser(description='Non-negative matrix factorization using squared Euclidean distance')					
	parser.add_argument('-kclusters', default='2', dest='kclusters', type=int, help='[INT] Number of subtypes or clusters to expect, must be smaller than m columns and n rows of input data')
	parser.add_argument('-input', required=True, dest='matrixFile', help='Full path to tab-delimited "matrix" file')
	parser.add_argument('-iterations', default='1000', dest='iterations', type=int, help='[INT] Number of iterations requried for convergence')
	parser.add_argument('--noPlotOut', default=True, dest='makePlot', action='store_false', help='[BOOLEAN] True or False, output image files of costfunction optimization')
	parser.add_argument('--colNames', default='noXLabels', dest='colNames', type=str, help='full path to file of sample names in order of matrix, one name per line')
	parser.add_argument('--rowNames', default='noYLabels', dest='rowNames', type=str, help='full path to file of feature/attribute names in order of matrix, one name per line')
	parser.add_argument('--output', default=os.getcwd(), dest='outPath', type=str, help='full path to output directory')
	args=parser.parse_args()
	
	uniqueName=strftime("%Y-%m-%d_%H:%M:%S", gmtime())
	#read in matrix file and randomly initialize matrix W and H
	observed, W, H=matrixInitialization(matrixFile=args.matrixFile, k=args.kclusters);
	
	#stores the quantification values of the cost function
	qualityApprox=[]
	iterConverge=[]
	predictionAccuracy=open(str(args.outPath)+'squared_euclidean_distance_at_each_iteration_k='+str(args.kclusters)+'_'+str(uniqueName)+'.txt', 'w')
	for x in range(0, args.iterations):
		W, H, predicted=updateRules(observed, W, H)
		squaredEucDist=costFunction(observed, predicted)
		qualityApprox.append(squaredEucDist)
		iterConverge.append(x)
		#predictionAccuracy.write(str(x)+'\t'(squaredEucDist)+'\n')
	
	matrix_visualization(W, H, yAxisNames=args.rowNames, xAxisNames=args.colNames)

	#output accuracy plots
	if args.makePlot==True:
		visualizeConvergenceAccuracy(qualityApprox, iterConverge);

	runInfo=open(str(args.outPath)+'run_metrics_k='+str(args.kclusters)+'_'+str(uniqueName)+'.txt', 'w')
	runInfo.write('number_of_iterations'+'\t'+str(args.iterations)+'\n')
	runInfo.write('number_of_clusters'+'\t'+str(args.kclusters)+'\n')
	runInfo.write('mean_squared_Euclidean_Distance'+'\t'+str(np.mean(qualityApprox))+'\n')
	runInfo.write('std_squared_Euclidean_Distance'+'\t'+str(np.std(qualityApprox))+'\n')
	runInfo.write('min_squared_Euclidean_Distance'+'\t'+str(np.min(qualityApprox))+'\n')
	runInfo.write('max_squared_Euclidean_Distance'+'\t'+str(np.max(qualityApprox))+'\n')	
	
	#outputs predicted W, H, and final predicted V, matrices in tab delimited format
	matrixH=open(str(args.outPath)+'matrixH_final_clusterXcolumn_k='+str(args.kclusters)+'_'+str(uniqueName)+'.txt', 'w')
	for x in range(0, len(H)):
		for z in range(0, len(H[0])-1):
			matrixH.write(str(H[x][z])+'\t')
		matrixH.write(str(H[x][len(H[0])-1])+'\n')

	matrixW=open(str(args.outPath)+'matrixW_final_rowXcluster_k='+str(args.kclusters)+'_'+str(uniqueName)+'.txt', 'w')
	for x in range(0, len(W)):
		for z in range(0, len(W[0])-1):
			matrixW.write(str(W[x][z])+'\t')
		matrixW.write(str(W[x][len(W[0])-1])+'\n')

	predictedMatrix=open(str(args.outPath)+'predicted_matrix_final_k='+str(args.kclusters)+'_'+str(uniqueName)+'.txt', 'w')
	for x in range(0, len(predicted)):
		for z in range(0, len(predicted[0])-1):
			predictedMatrix.write(str(predicted[x][z])+'\t')
		predictedMatrix.write(str(predicted[x][len(predicted[0])-1])+'\n')

