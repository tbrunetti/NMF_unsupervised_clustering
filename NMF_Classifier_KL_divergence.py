from time import gmtime, strftime
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os

def matrixInitialization(inputMatrix, k):
	#reads in matrix data
	observed=[]
	with open(inputMatrix) as input:
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
	predicted=np.dot(W, H)
	return observed, predicted, k, nrows, mcols, W, H

#updates matrix W by using multipicative update rule that minimizes KL divergence using nonincreasing rule below
def updateW(observed, predicted, k, nrows, mcols, W, H):
	for a in range(0, k):
		for i in range(0, nrows):
			numerator=0
			denominator=0
			for u in range(0, mcols):
				numerator=numerator+((H[a][u]*observed[i][u])/predicted[i][u])
				denominator=denominator+H[a][u]
			W[i][a]=W[i][a]*(numerator/denominator)
	predicted=np.dot(W, H)
	return W, predicted
	
#updates matrix H by using multipicative update rule that minimizes KL divergence using nonincreasing rule below
def updateH(observed, predicted, k, nrows, mcols, W, H):
	for a in range(0, k):
		for u in range(0, mcols):
			numerator=0
			denominator=0
			for i in range(0, nrows):
				numerator=numerator+((W[i][a]*observed[i][u])/predicted[i][u])
				denominator=denominator+W[i][a]
			H[a][u]=H[a][u]*(numerator/denominator)
	predicted=np.dot(W,H)
	return H, predicted

#measure of quality of approximation is determined by Kullback-Leibler divergence/relative entropy
def costFunction(observed, predicted):
	observedMultLogPredicted=np.multiply(observed, np.log(np.divide(observed, predicted)))
	#in case the log cannot be taken or index in predicted is 0 and cannot divide by zero,
	#replace all nan with 0
	observedMultLogPredicted=np.nan_to_num(observedMultLogPredicted)
	subtractObserved=np.subtract(observedMultLogPredicted, observed)
	addPredicted=np.add(subtractObserved, predicted)
	KL_divergence=np.sum(addPredicted)
	return KL_divergence

def visualizeConvergenceAccuracy(qualityApprox, iterConverge):
	#a line graph of the minimized costFunction across all iterations
	fig1, ax1=plt.subplots()
	plt.plot(iterConverge, qualityApprox)
	plt.title('Kullback-Leibler divergence between observed vs predicted')
	plt.axis([0, max(iterConverge), 0, max(qualityApprox)+50])   
	plt.ylabel('Kullback-Leibler divergence')
	plt.xlabel('iterations')
	plt.savefig(str(args.outPath)+'KLdiv_over_all_iterations_k='+str(args.kclusters)+'_'+str(uniqueName)+'.png')

	#zoomed in plot of the first 10 iterations
	plt.plot(iterConverge, qualityApprox)
	plt.title('Kullback-Leibler divergence between observed vs predicted')
	plt.axis([0, 10, 0, max(qualityApprox)+50])
	plt.ylabel('Kullback-Leibler divergence')
	plt.xlabel('iterations')
	plt.savefig(str(args.outPath)+'KLdiv_over_first_10_iterations_k='+str(args.kclusters)+'_'+str(uniqueName)+'.png')

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
			plt.savefig(str(args.outPath)+'KLdiv_matrixW_visualization_k='+str(args.kclusters)+'_'+str(uniqueName)+'.png')
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
			plt.savefig(str(args.outPath)+'KLdiv_matrixW_visualization_k='+str(args.kclusters)+'_'+str(uniqueName)+'.png')
	
	#construct heatmap fo matrix H
	def matrix_H():
	#if no names are provided for the columns it makes a heatmap with no labels
		if xAxisNames=='noXLabels':
			fig3, ax3=plt.subplots()
			heatmap=ax3.pcolor(H, cmap=plt.cm.RdYlGn)
			cbar = plt.colorbar(heatmap)
			plt.ylabel('gene expression profiles (k clusters)')
			plt.xlabel('sample ID')
			plt.savefig(str(args.outPath)+'KLdiv_matrixH_visualization_k='+str(args.kclusters)+'_'+str(uniqueName)+'.png')
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
			plt.savefig(str(args.outPath)+'KLdiv_matrixH_visualization_k='+str(args.kclusters)+'_'+str(uniqueName)+'.png')

	matrix_W();
	matrix_H();



if __name__=='__main__':
	#arguments and options for calling the classifier for KL-divergence and prediction of W, H
	parser=argparse.ArgumentParser(description='Non-negative matrix factorization using Kullback-Leibler divergence')
	parser.add_argument('-input', required=True, dest='matrixFile', help='Full path to tab-delimited "matrix" file')
	parser.add_argument('-kclusters', default='2', dest='kclusters', type=int, help='[INT] Number of subtypes or clusters to expect, must be smaller than m columns and n rows of input data')
	parser.add_argument('-iterations', default='1000', dest='iterations', type=int, help='[INT] Number of iterations requried for convergence')
	parser.add_argument('--noPlotOut', default=True, dest='makePlot', action='store_false', help='[BOOLEAN] True or False, output image files of costfunction optimization')
	parser.add_argument('--colNames', default='noXLabels', dest='colNames', type=str, help='full path to file of sample names in order of matrix, one name per line')
	parser.add_argument('--rowNames', default='noYLabels', dest='rowNames', type=str, help='full path to file of feature/attribute names in order of matrix, one name per line')
	parser.add_argument('--output', default=os.getcwd(), dest='outPath', type=str, help='full path to output directory')
	args=parser.parse_args()
	
	uniqueName=strftime("%Y-%m-%d_%H:%M:%S", gmtime())
	observed, predicted, k, nrows, mcols, W, H=matrixInitialization(inputMatrix=args.matrixFile, k=args.kclusters);
	
	qualityApprox=[]
	iterConverge=[]
	for i in range(0, args.iterations):
		W, predicted=updateW(observed, predicted, k, nrows, mcols, W, H);
		H, predicted=updateH(observed, predicted, k, nrows, mcols, W, H);
		KL_divergence=costFunction(observed, predicted)
		qualityApprox.append(KL_divergence)
		iterConverge.append(i)

	#calls to make visualizations
	if args.makePlot==True:
		visualizeConvergenceAccuracy(qualityApprox, iterConverge);
		matrix_visualization(W, H, yAxisNames=args.rowNames, xAxisNames=args.colNames)
	
	#creates a file with run statistics
	runInfo=open(str(args.outPath)+'KLdiv_run_metrics_k='+str(args.kclusters)+'_'+str(uniqueName)+'.txt', 'w')
	runInfo.write('number_of_iterations'+'\t'+str(args.iterations)+'\n')
	runInfo.write('number_of_clusters'+'\t'+str(args.kclusters)+'\n')
	runInfo.write('mean_Kullback-Leibler_divergence'+'\t'+str(np.mean(qualityApprox))+'\n')
	runInfo.write('std_Kullback-Leibler_divergence'+'\t'+str(np.std(qualityApprox))+'\n')
	runInfo.write('min_Kullback-Leibler_divergence'+'\t'+str(np.min(qualityApprox))+'\n')
	runInfo.write('max_Kullback-Leibler_divergence'+'\t'+str(np.max(qualityApprox))+'\n')	

	
	#outputs predicted W, H, and final predicted V, matrices in tab delimited format
	matrixH=open(str(args.outPath)+'KLdiv_matrixH_final_clusterXcolumn_k='+str(args.kclusters)+'_'+str(uniqueName)+'.txt', 'w')
	for x in range(0, len(H)):
		for z in range(0, len(H[0])-1):
			matrixH.write(str(H[x][z])+'\t')
		matrixH.write(str(H[x][len(H[0])-1])+'\n')

	matrixW=open(str(args.outPath)+'KLdiv_matrixW_final_rowXcluster_k='+str(args.kclusters)+'_'+str(uniqueName)+'.txt', 'w')
	for x in range(0, len(W)):
		for z in range(0, len(W[0])-1):
			matrixW.write(str(W[x][z])+'\t')
		matrixW.write(str(W[x][len(W[0])-1])+'\n')

	predictedMatrix=open(str(args.outPath)+'KLdiv_predicted_matrix_final_k='+str(args.kclusters)+'_'+str(uniqueName)+'.txt', 'w')
	for x in range(0, len(predicted)):
		for z in range(0, len(predicted[0])-1):
			predictedMatrix.write(str(predicted[x][z])+'\t')
		predictedMatrix.write(str(predicted[x][len(predicted[0])-1])+'\n')