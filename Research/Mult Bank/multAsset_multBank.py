
import numpy as np
# from numpy.linalg import inv
# from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from main_mult import *

#Input Initial Values Bank1
b1=.04	#constant in f(gamma(t))
s_initial1=np.array([[20.0],[11.0]])  #initial stocks
s1=np.copy(s_initial1) #array of current stock values (used in algorithm) 
alpha1=1.0
event1=-11
Gamma1= np.array([[0.0], [0.0]]) #Inital number of stocks sold at t. Always 0.
psi1= np.array([[0.0], [0.0]]) #Initial cash made from liquidation. Always 0.
ratio1=np.array([0]) #Start out selling 0% of assets

#Input Initial Values Bank 2
b2=.04	#constant in f(gamma(t))
s_initial2=np.array([[11.0],[11.0]])  #initial stocks
s2=np.copy(s_initial2) #array of current stock values (used in algorithm) 
alpha2=1.0
event2=100
Gamma2= np.array([[0.0], [0.0]]) #Inital number of stocks sold at t. Always 0.
psi2= np.array([[0.0], [0.0]]) #Initial cash made from liquidation. Always 0.
ratio2=np.array([0]) #Start out selling 0% of assets

#Store All Bank Values inside Arrays
allB= np.array([[b1], [b2]])
allS= np.array([s1, s2])
allS_init=np.array([s_initial1, s_initial2])
allStart= np.array([0, 0])
currentSeller= np.array([0, 0])
allAlpha= np.array([alpha1, alpha2])
allEvent= np.array([event1, event2])
allGamma= np.array([Gamma1, Gamma2])
allPsi= np.array([psi1, psi2])
allRatio= np.array([0, 0])

t=np.array([0.0])
dt=0.001
q = np.array([[1.0], [1.0]])  #Price of stock
theta_min= 0.08
k=0
assets=2 #Insert number of assets
banks=2

#Itialize F_t values
f_t= np.exp(-t[k])
plot_ft = np.array([np.exp(-t[k])])
fprime_t= -f_t

#Calculate the Cumulative Gamma values for each asset
allAssetsTotals= findAllGamma(allGamma,k, assets)

#Itialize F_g values
f_g=fGamma(allB, allAssetsTotals)
fprime_g=fGamma_prime(allB, f_g, assets)

#Initialize theta
theta_dot=updateThetas(allEvent, allPsi, allS_init, allRatio, q, allAlpha, banks, dt, allB)
allTheta=np.array(theta_dot)

while (t[k]<1): 
	print t[k]
	#Calculate the Cumulative Gamma values for each asset
	allAssetsTotals= findAllGamma(allGamma,k, assets)

	#update F values
	f_t= np.exp(-t[k])
	plot_ft = np.concatenate((plot_ft, [np.exp(-t[k])]))
	fprime_t= -f_t
	f_g=fGamma(allB, allAssetsTotals)
	fprime_g=fGamma_prime(allB, f_g, assets)

	#Calculate gamma, psi with previous k values for each bank 
	rowCounter=0
	for row in allGamma:
		Gamma_dot, Psi_dot = updateBank(allGamma[rowCounter], allStart[rowCounter], allS[rowCounter], allRatio[rowCounter], allPsi[rowCounter], q, dt)
		allGamma[rowCounter]=Gamma_dot
		allPsi[rowCounter]=Psi_dot
		if rowCounter==0: #UPDATE IF ADDING BANKS
			Gamma1= np.concatenate((Gamma1, Gamma_dot), axis=1)
			psi1=np.concatenate((psi1, Psi_dot), axis=1)
		else:
			Gamma2= np.concatenate((Gamma2, Gamma_dot), axis=1)
			psi2=np.concatenate((psi2, Psi_dot), axis=1)
		rowCounter = rowCounter+1

	#update q values
	q_dot=f_t*f_g
	q = np.concatenate((q, np.reshape(q_dot, (assets,1))), axis=1)

	#Calculate Theta last
	theta_dot=updateThetas(allEvent, allPsi, allS_init, allRatio, q, allAlpha, banks, dt, allTheta)
	allTheta=np.concatenate((allTheta, theta_dot), axis=1)

	#Update current stock units UPDATE IF ADDING ASSETS
	s_dot1= np.array(np.reshape(s1[:,0], (assets,1))-allGamma[0])
	s1= np.concatenate((s1, np.reshape(s_dot1, (assets,1))), axis=1)
	s_dot2= np.array(np.reshape(s2[:,0], (assets,1))-allGamma[1])
	s2= np.concatenate((s2, np.reshape(s_dot2, (assets,1))), axis=1)
	allS=np.array([s_dot1, s_dot2])

	#Update t
	t= np.concatenate((t, [t[k]+dt]))
	#Update k
	k=k+1
	
	m=0
	for row in allTheta:	
		if (allTheta[m][-1]-theta_min < 1e-20):
			np.put(allStart, m, 1)
		m=m+1

	minimum= np.amin(allTheta[:,k],0)

	#Check if below the threshold
	if (minimum-theta_min < 1e-20):
		#Set the start value to 1 to indicate selling, find ratio UPDATE IF ADDING BANKS
		m=0

		while m < len(allStart):
			if allStart[m]==1: #Update the ratio for selling banks  
				# newRatio=findRatio(allEvent[m], allPsi[m], q[:,k], theta_min, allAlpha[m], allS_init[m], k)
				newRatio=findRatioDeriv(event1, allPsi[m], q[:,k], allS_init[m], k, q_dot.item(0), allRatio[m], dt)
				allRatio[m]=newRatio
				if m==0:
					ratio1=np.concatenate((ratio1, [newRatio]))
				else:
					ratio2=np.concatenate((ratio2, [newRatio]))
			else:
				if m==0:
					ratio1=np.concatenate((ratio1, [ratio1.item((k-1))]))
				else:
					ratio2=np.concatenate((ratio2, [ratio2.item((k-1))]))
			m=m+1
		#Delete the theta, gamma, psi, and q values for the last event t UPDATE IF ADDING BANKS
		allTheta=np.delete(allTheta, k, 1)
		Gamma1=np.delete(Gamma1, k, 1)
		psi1=np.delete(psi1, k, 1)
		Gamma2=np.delete(Gamma2, k, 1)
		psi2=np.delete(psi2, k, 1)
		dummy1= np.array([[Gamma1[0][-1]], [Gamma1[1][-1]]])
		dummy2= np.array([[Gamma2[0][-1]], [Gamma2[1][-1]]])
		dummy3= np.array([[psi1[0][-1]], [psi1[1][-1]]])
		dummy4= np.array([[psi2[0][-1]], [psi2[1][-1]]])
		allGamma= np.array([dummy1, dummy2])
		allPsi= np.array([dummy3, dummy4])
		q=np.delete(q, k, 1)
		plot_ft=np.delete(plot_ft, k, 0)
		limit=t[k]
		t=np.delete(t, k, None)
		#Set the dt value to 100th of the original
		dt_new=dt/100

		#Update the itertion count to the previous
		k=k-1
		c=0


		while ((limit-t[k])>dt_new):  #iterate until "whole" interval 
			#Calculate the Cumulative Gamma values for each asset
			allAssetsTotals= findAllGamma(allGamma,k, assets)
			#update F values
			f_t= np.exp(-t[k])
			plot_ft = np.concatenate((plot_ft, [np.exp(-t[k])]))
			fprime_t= -f_t
			f_g=fGamma(allB, allAssetsTotals)
			fprime_g=fGamma_prime(allB, f_g, assets)

			#Calculate gamma, psi with previous k values for each bank 
			rowCounter=0

			for row in allGamma:
				Gamma_dot, Psi_dot = updateBank(allGamma[rowCounter], allStart[rowCounter], allS[rowCounter], allRatio[rowCounter], allPsi[rowCounter], q, dt)
				allGamma[rowCounter]=Gamma_dot
				allPsi[rowCounter]=Psi_dot
				if rowCounter==0: #UPDATE IF ADDING BANKS
					Gamma1= np.concatenate((Gamma1, Gamma_dot), axis=1)
					psi1=np.concatenate((psi1, Psi_dot), axis=1)
				else:
					Gamma2= np.concatenate((Gamma2, Gamma_dot), axis=1)
					psi2=np.concatenate((psi2, Psi_dot), axis=1)
				rowCounter = rowCounter+1

			#update q values
			q_dot=f_t*f_g
			q = np.concatenate((q, np.reshape(q_dot, (assets,1))), axis=1)

			#Calculate Theta last
			theta_dot=updateThetas(allEvent, allPsi, allS_init, allRatio, q, allAlpha, banks, dt_new, allTheta)
			allTheta=np.concatenate((allTheta, theta_dot), axis=1)

			#Update current stock units UPDATE IF ADDING BANKS
			s_dot1= np.array(np.reshape(s1[:,0], (assets,1))-allGamma[0])
			s1= np.concatenate((s1, np.reshape(s_dot1, (assets,1))), axis=1)
			s_dot2= np.array(np.reshape(s2[:,0], (assets,1))-allGamma[1])
			s2= np.concatenate((s2, np.reshape(s_dot2, (assets,1))), axis=1)
			allS=np.array([s_dot1, s_dot2])

			if c>0: #UPDATE IF ADDING BANKS
				ratio1=np.concatenate((ratio1, [ratio1.item((k))]))	
				ratio2=np.concatenate((ratio2, [ratio2.item((k))]))
				allRatio= np.array([ratio1.item((k)), ratio2.item((k))])
			t= np.concatenate((t, [t[k]+dt_new]))	
			c=c+1
			k=k+1

	else: #UPDATE IF ADDING BANKS
		ratio1=np.concatenate((ratio1, [ratio1.item((k-1))]))
		ratio2=np.concatenate((ratio2, [ratio2.item((k-1))]))
		allRatio= np.array([ratio1.item((k-1)), ratio2.item((k-1))])


#Print final values and graphs
print("Final Theta Bank 1:")
print (allTheta[0][-1])
print("Final Theta Bank 2:")
print (allTheta[1][-1])
print "Ratio Bank 1"
print ratio1[-1]
print "Ratio Bank 2"
print ratio2[-1]
print("Final t:")
print(t[k])
print("Final S bank 1:")
print (allS[0])
print("Final S bank 2:")
print(allS[1])
tnew=np.reshape(t, ((k+1),1))
assetDummy=0
trows= len(t)
print("Final Gamma Bank 1 :")
print(Gamma1[0][-1])
print(Gamma1[1][-1])
print("Final Gamma Bank 2 :")
print(Gamma2[0][-1])
print(Gamma2[1][-1])
print("Final Psi Bank 1:")
print(psi1[0][-1])
print(psi1[1][-1])
print("Final Psi Bank 2:")
print(psi2[0][-1])
print(psi2[1][-1])
while (assetDummy < assets):
	print("Final q Asset " + repr(assetDummy + 1) + ":")
	print(q.item((assetDummy,k)))
	assetDummy=assetDummy+1
plt.figure(1)
plt.plot(tnew, allTheta[0], color='orange') 
# plt.plot(tnew, allTheta[1], color='blue') 
plt.plot(t,-ratio1, color='red')  
# plt.plot(t,-ratio2, color='green')  		
counter=1
for row in q:
	if counter==1:
		plt.plot(t, row, color='plum')
	else:
		plt.plot(t, row, color='purple')
	counter=counter+1
plt.title("Single-Bank, Multi-Asset Method 2 Results")
plt.xlabel("Time (t)")
plt.ylabel("Scaled Simulation Value")
plum1_patch = mpatches.Patch(color='plum', label='Asset 1 Price')
plum2_patch = mpatches.Patch(color='purple', label='Asset 2 Price')
orange_patch = mpatches.Patch(color='orange', label='Capital Ratio')
# blue_patch = mpatches.Patch(color='blue', label='Capital Ratio Bank 2')
red_patch = mpatches.Patch(color='red', label='Selling Ratio')
# green_patch = mpatches.Patch(color='green', label='Selling Ratio Bank 2')
plt.legend(handles=[plum1_patch, plum2_patch, orange_patch, red_patch])
plt.show()
plt.figure(2)
plt.plot(tnew, allTheta[0], color='orange') 
plt.plot(tnew, allTheta[1], color='blue') 
plt.plot(t,-ratio1, color='red')  
plt.plot(t,-ratio2, color='green')		
plt.title("Single-Bank, Multi-Asset Method 2 Selling and Capital Ratio Results")
plt.xlabel("Time (t)")
plt.ylabel("Scaled Simulation Value")
plt.legend(handles=[orange_patch, red_patch])
plt.show()
plt.figure(3)  		
counter=1
for row in q:
	if counter==1:
		plt.plot(t, row, color='red')
	else:
		plt.plot(t, row, color='blue')
	counter=counter+1
plt.plot(t, plot_ft, color='green')
plt.title("Price Impact on Asset Value")
plt.xlabel("Time (t)")
plt.ylabel("Scaled Simulation Value")
red_patch = mpatches.Patch(color='red', label='Asset 1 Price')
green_patch = mpatches.Patch(color='green', label='Asset Decay Without Price Impact')
blue_patch = mpatches.Patch(color='blue', label='Asset 2 Price')

plt.legend(handles=[red_patch, blue_patch, green_patch])
plt.show()