import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

def matrixInitialization():
	#reads in matrix data
	observed=[]
	with open(sys.argv[1]) as input:
		for line in input:
			line=line.rstrip('\n')
			line=line.split('\t')
			converted=[float(x) for x in line]
			observed.append(converted)

	observed=np.array(observed)
	nrows=len(observed)
	mcols=len(observed[0])

	#randomly initialize matrix W and H to small non-negative numbers
	W=np.random.rand(nrows, k)
	H=np.random.rand(k, mcols)
	predicted=np.dot(W, H)
	return observed, predicted, k, nrows, mcols, W, H

def updateW(observed, predicted, k, nrows, mcols, W, H):
	for a in range(0, k):
		for i in range(0, nrows):
			numerator=0
			denominator=0
			for u in range(0, mcols):
				numerator=numerator+((H[a][u]*V[i][u])/predicted[i][u])
				denominator=denominator+H[a][u]
			W[i][a]=W[i][a]*(numerator/denominator)
	predicted=np.dot(W, H)
	return W, predicted
	
def updateH(observed, predicted, k, nrows, mcols, W, H):
	for a in range(0, k):
		for u in range(0, mcols):
			numerator=0
			denominator=0
			for i in range(0, nrows):
				numerator=numerator+((W[i][a]*V[i][u])/predicted[i][u])
				denominator=denominator+W[i][a]
			H[a][u]=H[a][u]*(numerator/denominator)
	predicted=np.dot(W,H)
	return H, predicted

def costFunction():
	pass;

if __name__=='__main__':
	observed, predicted, k, nrows, mcols, W, H=matrixInitialization();
