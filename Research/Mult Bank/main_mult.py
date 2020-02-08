
import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA


#Functions to be called from the multi bank model scripts

def updateGamma(start, s, ratio):
	if (start==1):
		Gamma_dot=-ratio*s 
	else:
		Gamma_dot= 0
	return Gamma_dot

def findRatio(event, psi, q, theta_min, alpha, s, k):
	psiSum=np.cumsum(psi, axis=0)
	i=0
	for row in psi:
		if i==0:
			sqArray=np.array([s[i]*q[i]*(theta_min*alpha-1)])
		else:
			sqArray=np.concatenate((sqArray,[s[i]*q[i]*(theta_min*alpha-1)]))
		i=i+1
	psiSum=psiSum.item(((i-1),0))
	numerator=event+psiSum
	sqSum=np.cumsum(sqArray, axis=0)
	denominator=sqSum[i-1]
	ratio=1-(numerator/denominator)
	ratio=ratio[0]
	return ratio

def diversification(banks, allDegree):
	counter = 0
	sumD = 0
	while counter < banks:
		sumD = sumD + allDegree[counter]
		counter = counter+1
	return ((1/float(banks))*sumD)

def assetTotal(assets, allS_init):
	counter = 0
	asset1_total = 0
	asset2_total = 0
	asset3_total = 0
	asset4_total = 0
	asset5_total = 0
	asset6_total = 0
	#update if more assets
	while counter < assets :
		i = 0
		for row in allS_init:
			if allS_init[i][counter] != 0:
				if counter == 0:
					asset1_total = asset1_total + 1
				if counter == 1:
					asset2_total = asset2_total + 1
				if counter == 2:
					asset3_total = asset3_total + 1
				if counter == 3:
					asset4_total = asset4_total + 1
				if counter == 4:
					asset5_total = asset5_total + 1
				if counter == 5:
					asset6_total = asset6_total + 1
			i=i+1
		counter = counter + 1

	return np.array([asset1_total,asset2_total, asset3_total,asset4_total,asset5_total, asset6_total])

def assetDiversification(assets, assetsInPortfolios):
	Sum= np.cumsum(assetsInPortfolios, axis=0)
	return (1/float(assets)*Sum[-1])
	
def findAllGamma(allGamma, k, assets):
	allAssets = {}
	bankCounter=1
	for row in allGamma: #iterate through banks gamma arrays
		counter=0
		while counter < assets: #iterate through each banks gamma values (for each asset)
			allAssets[str(bankCounter)+"_"+str(counter+1)]=row[counter] # Key = assset #; Value = Gamma
			counter=counter+1
		bankCounter=bankCounter+1

	asset1_Total=0
	asset2_Total=0
	asset3_Total=0
	asset4_Total=0
	asset5_Total=0
	asset6_Total=0
	#UPDATE IF ASSETS ADDED
	for x in allAssets:
		if "_1" in x:
			asset1_Total=asset1_Total+allAssets[x]
		if "_2" in x:
			asset2_Total=asset2_Total+allAssets[x]	
		if "_3" in x:
			asset3_Total=asset3_Total+allAssets[x]	
		if "_4" in x:
			asset4_Total=asset4_Total+allAssets[x]
		if "_5" in x:
			asset5_Total=asset5_Total+allAssets[x]	
		if "_6" in x:
			asset6_Total=asset6_Total+allAssets[x]
	allAssetsTotals= np.array([[asset1_Total], [asset2_Total],[asset3_Total],[asset4_Total], [asset5_Total],[asset6_Total]])
	return allAssetsTotals

def fGamma(b, Gamma): 
	counter=0
	for row in Gamma:
		if counter == 0 :
			f_g=np.array([np.exp(-b[counter]*Gamma[counter])])
		else:
			f_g=np.concatenate((f_g, [np.exp(-b[counter]*Gamma[counter])]))
		counter=counter+1
	return f_g

def fGamma_prime(b, f_g, assets):
	j=0
	while (j<assets):
		if (j==0):
			fprime_g=np.array([-b[j]*f_g[j]])
		else:
			fprime_g=np.concatenate((fprime_g,[-b[j]*f_g[j]]))
		j=j+1
	return fprime_g

def updateTheta(event, psi, s, ratio, q, alpha, dt, allTheta):
	n=0
	for row in psi:
		sum1=[psi[n]+s[n]*(1-ratio)*q[n][-1]]
		sum2=[alpha*(s[n]*(1-ratio)*q[n][-1])]
		if sum1<0:
			sum1=[0]
		if n==0:
			theta1=np.array(sum1)
			theta2=np.array(sum2)
		else:
			theta1=np.concatenate((theta1, sum1))
			theta2=np.concatenate((theta2, sum2))
		n=n+1
	numSum= np.cumsum(theta1, axis=0)
	denSum= np.cumsum(theta2, axis=0)
	numerator= event+numSum[n-1]
	if numerator <0 :
		numerator=0
	denominator=denSum[n-1]
	theta_dot=numerator/denominator
	theta_dot=np.reshape(theta_dot, (1,1))
	return theta_dot

def updateThetas(allEvent, allPsi, allS, allRatio, q, allAlpha, banks, dt, allTheta):
	bankCounter=0
	while bankCounter < banks:
		if bankCounter==0:
			theta_dot=np.array(updateTheta(allEvent[bankCounter], allPsi[bankCounter], allS[bankCounter], allRatio[bankCounter], q, allAlpha[bankCounter], dt, allTheta[bankCounter]))
		else:
			theta_dot=np.concatenate((theta_dot,updateTheta(allEvent[bankCounter], allPsi[bankCounter], allS[bankCounter], allRatio[bankCounter], q, allAlpha[bankCounter], dt, allTheta[bankCounter])))
		bankCounter=bankCounter+1
	return theta_dot


def updateBank(Gamma, start, s, ratio, psi, q, dt):
	n=0 
	for row in Gamma:
		Gamma_new=updateGamma(start, s.item(n), ratio)
		if n==0:
			Gamma_dot=np.array([Gamma[n]+Gamma_new*dt])
			Psi_dot=np.array([psi[n]+Gamma_new*q[n][-1]*dt])
		else:
			Gamma_dot=np.concatenate((Gamma_dot, [Gamma[n]+Gamma_new*dt]))
			Psi_dot=np.concatenate((Psi_dot, [psi[n]+Gamma_new*q[n][-1]*dt]))

		n=n+1
	#Reshape arrays and append new values
	Gamma_dot=(np.reshape(Gamma_dot, (n,1))) 
	Psi_dot=(np.reshape(Psi_dot, (n,1)))
	return Gamma_dot, Psi_dot


