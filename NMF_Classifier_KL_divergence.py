import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

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

if __name__=='__main__':
	#arguments and options for calling the classifier for KL-divergence and prediction of W, H
	parser=argparse.ArgumentParser(description='Non-negative matrix factorization using Kullback-Leibler divergence')
	parser.add_argument('-input', required=True, dest='matrixFile', help='Full path to tab-delimited "matrix" file')
	parser.add_argument('-kclusters', default='2', dest='kclusters', type=int, help='[INT] Number of subtypes or clusters to expect, must be smaller than m columns and n rows of input data')
	parser.add_argument('-iterations', default='1000', dest='iterations', type=int, help='[INT] Number of iterations requried for convergence')
	parser.add_argument('--noPlotOut', default=True, dest='makePlot', action='store_false', help='[BOOLEAN] True or False, output image files of costfunction optimization')
	parser.add_argument('--colNames', default='noXLabels', dest='colNames', type=str, help='full path to file of sample names in order of matrix, one name per line')
	parser.add_argument('--rowNames', default='noYLabels', dest='rowNames', type=str, help='full path to file of feature/attribute names in order of matrix, one name per line')
	args=parser.parse_args()
	

	observed, predicted, k, nrows, mcols, W, H=matrixInitialization(inputMatrix=args.matrixFile, k=args.kclusters);
	for i in range(0, args.iterations):
		W, predicted=updateW(observed, predicted, k, nrows, mcols, W, H);
		H, predicted=updateH(observed, predicted, k, nrows, mcols, W, H);
		KL_divegence=costFunction(observed, predicted)
		print KL_divegence