
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from main_mult import *


#Multibank model implemented with Greece bank figures. Ran as a Use Case.

#Input Initial Values Bank1
b1=.04	#constant in f(gamma(t))
s_initial1=np.array([[12*.03],[17*.51],[13*.05],[10*.13],[3*.11],[1.5*.6]])  #initial stocks
s1=np.copy(s_initial1) #array of current stock values (used in algorithm) 
alpha1=1.0
event1=-5
Gamma1= np.array([[0.0], [0.0], [0.0],[0.0], [0.0], [0.0]]) #Inital number of stocks sold at t. Always 0.
psi1= np.array([[0.0], [0.0], [0.0],[0.0], [0.0], [0.0]]) #Initial cash made from liquidation. Always 0.
ratio1=np.array([0]) #Start out selling 0% of assets

#Input Initial Values Bank 2
b2=.04	#constant in f(gamma(t))
s_initial2=np.array([[12*.17],[17*.01],[0.0],[10*0.0],[3*.66],[1.5*.04]])  #initial stocks
s2=np.copy(s_initial2) #array of current stock values (used in algorithm) 
alpha2=1.0
event2=-10
Gamma2= np.array([[0.0], [0.0], [0.0],[0.0], [0.0], [0.0]]) #Inital number of stocks sold at t. Always 0.
psi2= np.array([[0.0], [0.0], [0.0],[0.0], [0.0], [0.0]]) #Initial cash made from liquidation. Always 0.
ratio2=np.array([0]) #Start out selling 0% of assets

#Input Initial Values Bank 3
b3=.04	#constant in f(gamma(t))
s_initial3=np.array([[12*0.0],[17*0.0],[0.0],[10*.83],[3*.01],[1.5*.03]])  #initial stocks
s3=np.copy(s_initial3) #array of current stock values (used in algorithm) 
alpha3=1.0
event3=-6
Gamma3= np.array([[0.0], [0.0], [0.0],[0.0], [0.0], [0.0]]) #Inital number of stocks sold at t. Always 0.
psi3= np.array([[0.0], [0.0], [0.0],[0.0], [0.0], [0.0]]) #Initial cash made from liquidation. Always 0.
ratio3=np.array([0]) #Start out selling 0% of assets

#Input Initial Values Bank 4
b4=.04	#constant in f(gamma(t))
s_initial4=np.array([[12*.05],[17*0.0],[11.96],[0.0],[0.0],[1.5*.17]])  #initial stocks
s4=np.copy(s_initial4) #array of current stock values (used in algorithm) 
alpha4=1.0
event4=-6
Gamma4= np.array([[0.0], [0.0], [0.0],[0.0], [0.0], [0.0]]) #Inital number of stocks sold at t. Always 0.
psi4= np.array([[0.0], [0.0], [0.0],[0.0], [0.0], [0.0]]) #Initial cash made from liquidation. Always 0.
ratio4=np.array([0]) #Start out selling 0% of assets

#Input Initial Values Bank 5
b5=.04	#constant in f(gamma(t))
s_initial5=np.array([[12*.01],[17*.15],[0.0],[0.4],[3*.03],[1.5*.02]])  #initial stocks
s5=np.copy(s_initial5) #array of current stock values (used in algorithm) 
alpha5=1.0
event5=-7
Gamma5= np.array([[0.0], [0.0], [0.0],[0.0], [0.0], [0.0]]) #Inital number of stocks sold at t. Always 0.
psi5= np.array([[0.0], [0.0], [0.0],[0.0], [0.0], [0.0]]) #Initial cash made from liquidation. Always 0.
ratio5=np.array([0]) #Start out selling 0% of assets

#Input Initial Values Bank 6
b6=.04	#constant in f(gamma(t))
s_initial6=np.array([[12*.72],[17*.32],[0.13],[0.0],[3*.16],[1.5*.03]])  #initial stocks
s6=np.copy(s_initial6) #array of current stock values (used in algorithm) 
alpha6=1.0
event6=-4
Gamma6= np.array([[0.0], [0.0], [0.0],[0.0], [0.0], [0.0]]) #Inital number of stocks sold at t. Always 0.
psi6= np.array([[0.0], [0.0], [0.0],[0.0], [0.0], [0.0]]) #Initial cash made from liquidation. Always 0.
ratio6=np.array([0]) #Start out selling 0% of assets

#Store All Bank Values inside Arrays
allB= np.array([[b1], [b2],[b3], [b4],[b5], [b6]])
allS= np.array([s1, s2, s3, s4, s5, s6])
allS_init=np.array([s_initial1, s_initial2, s_initial3, s_initial4,s_initial5, s_initial6])
allStart= np.array([0, 0,0,0,0,0])
allAlpha= np.array([alpha1, alpha2, alpha3, alpha4, alpha5, alpha6])
allEvent= np.array([event1, event2, event3, event4, event5, event6])
allGamma= np.array([Gamma1, Gamma2, Gamma3, Gamma4, Gamma5, Gamma6])
allPsi= np.array([psi1, psi2, psi3, psi4, psi5, psi6])
allRatio= np.array([0, 0, 0, 0, 0, 0])
allDegree= np.array([6, 6, 6, 6, 6, 5])
t=np.array([0.0])
plot_t=np.array([0.0])
dt=0.01
theta_min= 0.08
k=0
assets=6
banks=6

System_Diversification_Banks = diversification(banks, allDegree) 
assetsInPortfolios = assetTotal(assets, allS_init)
System_Diversification_Assets= assetDiversification(assets, assetsInPortfolios)
Bank_Diversification = allDegree / float(assets)
Asset_Diversification = assetsInPortfolios / float(banks)



#Itialize F_t values
f_t= np.exp(-t[k])
plot_ft = np.array([np.exp(-t[k])])
fprime_t= -f_t

# #Calculate the Cumulative Gamma values for each asset
allAssetsTotals= findAllGamma(allGamma,k, assets)
assetTotals = assetTotal(assets, allS_init)

#Itialize F_g values
f_g=fGamma(allB, allAssetsTotals)
fprime_g=fGamma_prime(allB, f_g, assets)

q_dot=f_t*f_g
q = np.array(np.reshape(q_dot, (assets,1)))

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
		if rowCounter==0:
			Gamma1= np.concatenate((Gamma1, Gamma_dot), axis=1)
			psi1=np.concatenate((psi1, Psi_dot), axis=1)
		if rowCounter==1:
			Gamma2= np.concatenate((Gamma2, Gamma_dot), axis=1)
			psi2=np.concatenate((psi2, Psi_dot), axis=1)
		if rowCounter==2:
			Gamma3= np.concatenate((Gamma3, Gamma_dot), axis=1)
			psi3=np.concatenate((psi3, Psi_dot), axis=1)
		if rowCounter==3:
			Gamma4= np.concatenate((Gamma4, Gamma_dot), axis=1)
			psi4=np.concatenate((psi4, Psi_dot), axis=1)
		if rowCounter==4:
			Gamma5= np.concatenate((Gamma5, Gamma_dot), axis=1)
			psi5=np.concatenate((psi5, Psi_dot), axis=1)
		else:
			Gamma6= np.concatenate((Gamma6, Gamma_dot), axis=1)
			psi6=np.concatenate((psi6, Psi_dot), axis=1)
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
	s_dot3= np.array(np.reshape(s3[:,0], (assets,1))-allGamma[2])
	s3= np.concatenate((s3, np.reshape(s_dot3, (assets,1))), axis=1)
	s_dot4= np.array(np.reshape(s4[:,0], (assets,1))-allGamma[3])
	s4= np.concatenate((s4, np.reshape(s_dot4, (assets,1))), axis=1)
	s_dot5= np.array(np.reshape(s5[:,0], (assets,1))-allGamma[4])
	s5= np.concatenate((s5, np.reshape(s_dot5, (assets,1))), axis=1)
	s_dot6= np.array(np.reshape(s6[:,0], (assets,1))-allGamma[5])
	s6= np.concatenate((s6, np.reshape(s_dot6, (assets,1))), axis=1)
	allS=np.array([s_dot1, s_dot2, s_dot3,s_dot4, s_dot5, s_dot6])

	t= np.concatenate((t, [t[k]+dt]))
	plot_t= np.concatenate((plot_t, [plot_t[-1]+dt]))
	k=k+1
	
	m=0
	for row in allTheta:	
		if (allTheta[m][-1]-theta_min < 1e-20):
			np.put(allStart, m, 1)
		m=m+1

	minimum= np.amin(allTheta[:,k],0)

	#Check if below the threshold
	if (minimum-theta_min < 1e-20):
		#Set the start value to 1 to indicate selling, find ratio
		m=0
		while m < len(allStart):
			if allStart[m]==1: #Update the ratio for selling banks  
				newRatio=findRatio(allEvent[m], allPsi[m], q[:,k], theta_min, allAlpha[m], allS_init[m], k)
				allRatio[m]=newRatio
				if m==0:
					ratio1=np.concatenate((ratio1, [newRatio]))
				if m==1:
					ratio2=np.concatenate((ratio2, [newRatio]))
				if m==2:
					ratio3=np.concatenate((ratio3, [newRatio]))
				if m==3:
					ratio4=np.concatenate((ratio4, [newRatio]))
				if m==4:
					ratio5=np.concatenate((ratio5, [newRatio]))
				if m==5:
					ratio6=np.concatenate((ratio6, [newRatio]))
			else:
				if m==0:
					ratio1=np.concatenate((ratio1, [ratio1.item((k-1))]))
				if m==1:
					ratio2=np.concatenate((ratio2, [ratio2.item((k-1))]))
				if m==2:
					ratio3=np.concatenate((ratio3, [ratio3.item((k-1))]))
				if m==3:
					ratio4=np.concatenate((ratio4, [ratio4.item((k-1))]))
				if m==4:
					ratio5=np.concatenate((ratio5, [ratio5.item((k-1))]))
				if m==5:
					ratio6=np.concatenate((ratio6, [ratio6.item((k-1))]))
			m=m+1

		#Delete the theta, gamma, psi, and q values for the last event t
		allTheta=np.delete(allTheta, k, 1)
		Gamma1=np.delete(Gamma1, k, 1)
		psi1=np.delete(psi1, k, 1)
		Gamma2=np.delete(Gamma2, k, 1)
		psi2=np.delete(psi2, k, 1)

		#Update if adding assets
		dummy1= np.array([[Gamma1[0][-1]], [Gamma1[1][-1]], [Gamma1[2][-1]],[Gamma1[3][-1]], [Gamma1[4][-1]], [Gamma1[5][-1]]])
		dummy2= np.array([[Gamma2[0][-1]], [Gamma2[1][-1]], [Gamma2[2][-1]],[Gamma2[3][-1]], [Gamma2[4][-1]], [Gamma2[5][-1]]])
		dummy3= np.array([[psi1[0][-1]], [psi1[1][-1]], [psi1[2][-1]],[psi1[3][-1]], [psi1[4][-1]], [psi1[5][-1]]])
		dummy4= np.array([[psi2[0][-1]], [psi2[1][-1]], [psi2[2][-1]],[psi2[3][-1]], [psi2[4][-1]], [psi2[5][-1]]])

		dummy5= np.array([[Gamma3[0][-1]], [Gamma3[1][-1]], [Gamma3[2][-1]],[Gamma3[3][-1]], [Gamma3[4][-1]], [Gamma3[5][-1]]])
		dummy6= np.array([[Gamma4[0][-1]], [Gamma4[1][-1]], [Gamma4[2][-1]],[Gamma4[3][-1]], [Gamma4[4][-1]], [Gamma4[5][-1]]])
		dummy7= np.array([[psi3[0][-1]], [psi3[1][-1]], [psi3[2][-1]],[psi3[3][-1]], [psi3[4][-1]], [psi3[5][-1]]])
		dummy8= np.array([[psi4[0][-1]], [psi4[1][-1]], [psi4[2][-1]],[psi4[3][-1]], [psi4[4][-1]], [psi4[5][-1]]])

		dummy9= np.array([[Gamma5[0][-1]], [Gamma5[1][-1]], [Gamma5[2][-1]],[Gamma5[3][-1]], [Gamma5[4][-1]], [Gamma5[5][-1]]])
		dummy10= np.array([[Gamma6[0][-1]], [Gamma6[1][-1]], [Gamma6[2][-1]],[Gamma6[3][-1]], [Gamma6[4][-1]], [Gamma6[5][-1]]])
		dummy11= np.array([[psi5[0][-1]], [psi5[1][-1]], [psi5[2][-1]],[psi5[3][-1]], [psi5[4][-1]], [psi5[5][-1]]])
		dummy12= np.array([[psi6[0][-1]], [psi6[1][-1]], [psi6[2][-1]],[psi6[3][-1]], [psi6[4][-1]], [psi6[5][-1]]])

		allGamma= np.array([dummy1, dummy2, dummy5, dummy6, dummy9, dummy10])
		allPsi= np.array([dummy3, dummy4, dummy7, dummy8, dummy11, dummy12])
		q=np.delete(q, k, 1)
		# plot_ft=np.delete(plot_ft, k, 0)
		limit=t[k]
		t=np.delete(t, k, None)

		dt_new=dt/100
		k=k-1
		c=0


		while ((limit-t[k])>dt_new):  #iterate until "whole" interval 
			allAssetsTotals= findAllGamma(allGamma,k, assets)

			f_t= np.exp(-t[k])
			# plot_ft = np.concatenate((plot_ft, [np.exp(-t[k])]))
			fprime_t= -f_t
			f_g=fGamma(allB, allAssetsTotals)
			fprime_g=fGamma_prime(allB, f_g, assets)

			#Calculate gamma, psi with previous k values for each bank 
			rowCounter=0
			for row in allGamma:
				Gamma_dot, Psi_dot = updateBank(allGamma[rowCounter], allStart[rowCounter], allS[rowCounter], allRatio[rowCounter], allPsi[rowCounter], q, dt)
				allGamma[rowCounter]=Gamma_dot
				allPsi[rowCounter]=Psi_dot
				if rowCounter==0:
					Gamma1= np.concatenate((Gamma1, Gamma_dot), axis=1)
					psi1=np.concatenate((psi1, Psi_dot), axis=1)
				if rowCounter==1:
					Gamma2= np.concatenate((Gamma2, Gamma_dot), axis=1)
					psi2=np.concatenate((psi2, Psi_dot), axis=1)
				if rowCounter==2:
					Gamma3= np.concatenate((Gamma3, Gamma_dot), axis=1)
					psi3=np.concatenate((psi3, Psi_dot), axis=1)
				if rowCounter==3:
					Gamma4= np.concatenate((Gamma4, Gamma_dot), axis=1)
					psi4=np.concatenate((psi4, Psi_dot), axis=1)
				if rowCounter==4:
					Gamma5= np.concatenate((Gamma5, Gamma_dot), axis=1)
					psi5=np.concatenate((psi5, Psi_dot), axis=1)
				if rowCounter==5:
					Gamma6= np.concatenate((Gamma6, Gamma_dot), axis=1)
					psi6=np.concatenate((psi6, Psi_dot), axis=1)
				rowCounter = rowCounter+1

			#update q values
			q_dot=f_t*f_g
			q = np.concatenate((q, np.reshape(q_dot, (assets,1))), axis=1)

			#Calculate Theta last
			theta_dot=updateThetas(allEvent, allPsi, allS_init, allRatio, q, allAlpha, banks, dt_new, allTheta)
			allTheta=np.concatenate((allTheta, theta_dot), axis=1)

			#Update current stock units UPDATE IF ADDING ASSETS
			s_dot1= np.array(np.reshape(s1[:,0], (assets,1))-allGamma[0])
			s1= np.concatenate((s1, np.reshape(s_dot1, (assets,1))), axis=1)
			s_dot2= np.array(np.reshape(s2[:,0], (assets,1))-allGamma[1])
			s2= np.concatenate((s2, np.reshape(s_dot2, (assets,1))), axis=1)
			s_dot3= np.array(np.reshape(s3[:,0], (assets,1))-allGamma[2])
			s3= np.concatenate((s3, np.reshape(s_dot3, (assets,1))), axis=1)
			s_dot4= np.array(np.reshape(s4[:,0], (assets,1))-allGamma[3])
			s4= np.concatenate((s4, np.reshape(s_dot4, (assets,1))), axis=1)
			s_dot5= np.array(np.reshape(s5[:,0], (assets,1))-allGamma[4])
			s5= np.concatenate((s5, np.reshape(s_dot5, (assets,1))), axis=1)
			s_dot6= np.array(np.reshape(s6[:,0], (assets,1))-allGamma[5])
			s6= np.concatenate((s6, np.reshape(s_dot6, (assets,1))), axis=1)
			allS=np.array([s_dot1, s_dot2, s_dot3,s_dot4, s_dot5, s_dot6])

			if c>0:
				ratio1=np.concatenate((ratio1, [ratio1.item((k))]))	
				ratio2=np.concatenate((ratio2, [ratio2.item((k))]))
				ratio3=np.concatenate((ratio3, [ratio3.item((k))]))	
				ratio4=np.concatenate((ratio4, [ratio4.item((k))]))
				ratio5=np.concatenate((ratio5, [ratio5.item((k))]))	
				ratio6=np.concatenate((ratio6, [ratio6.item((k))]))
				allRatio= np.array([ratio1.item((k)), ratio2.item((k)), ratio3.item((k)), ratio4.item((k)), ratio5.item((k)), ratio6.item((k))])
			t= np.concatenate((t, [t[k]+dt_new]))	
			c=c+1
			k=k+1

	else:
		ratio1=np.concatenate((ratio1, [ratio1.item((k-1))]))
		ratio2=np.concatenate((ratio2, [ratio2.item((k-1))]))
		ratio3=np.concatenate((ratio3, [ratio3.item((k-1))]))
		ratio4=np.concatenate((ratio4, [ratio4.item((k-1))]))
		ratio5=np.concatenate((ratio5, [ratio5.item((k-1))]))
		ratio6=np.concatenate((ratio6, [ratio6.item((k-1))]))
		allRatio= np.array([ratio1.item((k-1)), ratio2.item((k-1)), ratio3.item((k-1)), ratio4.item((k-1)), ratio5.item((k-1)), ratio6.item((k-1))])









#PRINTING AND PLOTTING
print("Final Theta Bank 1:")
print (allTheta[0][-1])
print("Final Theta Bank 2:")
print (allTheta[1][-1])
print("Final Theta Bank 3:")
print (allTheta[2][-1])
print("Final Theta Bank 4:")
print (allTheta[3][-1])
print("Final Theta Bank 5:")
print (allTheta[4][-1])
print("Final Theta Bank 6:")
print (allTheta[5][-1])
print "Ratio Bank 1"
print -ratio1[-1]
print "Ratio Bank 2"
print -ratio2[-1]
print "Ratio Bank 3"
print -ratio3[-1]
print "Ratio Bank 4"
print -ratio4[-1]
print "Ratio Bank 5"
print -ratio5[-1]
print "Ratio Bank 6"
print -ratio6[-1]
print("Final t:")
print(t[k])
print("Final S bank 1:")
print (allS[0])
print("Final S bank 2:")
print(allS[1])
print("Final S bank 3:")
print (allS[2])
print("Final S bank 4:")
print(allS[3])
print("Final S bank 5:")
print (allS[4])
print("Final S bank 6:")
print(allS[5])

print "--System Performance Review-- "
print "Number of Failing Firms:"
print np.count_nonzero(allStart == 1)
print "System_Diversification_Banks:" 
print System_Diversification_Banks 
print "assetsInPortfolios:"
print assetsInPortfolios 
print "System_Diversification_Assets:"
print System_Diversification_Assets
print"Bank_Diversification:"
print Bank_Diversification 
print "Asset_Diversification:"
print Asset_Diversification 

print("Final Gamma Bank 1 :")
print(Gamma1[0][-1])
print(Gamma1[1][-1])
print(Gamma1[2][-1])
print(Gamma1[3][-1])
print(Gamma1[4][-1])
print(Gamma1[5][-1])
print("Final Gamma Bank 2 :")
print(Gamma2[0][-1])
print(Gamma2[1][-1])
print(Gamma2[2][-1])
print(Gamma2[3][-1])
print(Gamma2[4][-1])
print(Gamma2[5][-1])
print("Final Gamma Bank 3 :")
print(Gamma3[0][-1])
print(Gamma3[1][-1])
print(Gamma3[2][-1])
print(Gamma3[3][-1])
print(Gamma3[4][-1])
print(Gamma3[5][-1])
print("Final Gamma Bank 4 :")
print(Gamma4[0][-1])
print(Gamma4[1][-1])
print(Gamma4[2][-1])
print(Gamma4[3][-1])
print(Gamma4[4][-1])
print(Gamma4[5][-1])
print("Final Gamma Bank 5 :")
print(Gamma5[0][-1])
print(Gamma5[1][-1])
print(Gamma5[2][-1])
print(Gamma5[3][-1])
print(Gamma5[4][-1])
print(Gamma5[5][-1])
print("Final Gamma Bank 6 :")
print(Gamma6[0][-1])
print(Gamma6[1][-1])
print(Gamma6[2][-1])
print(Gamma6[3][-1])
print(Gamma6[4][-1])
print(Gamma6[5][-1])

# print("Final Psi Bank 1:")
# print(psi1[0][-1])
# print(psi1[1][-1])
# print("Final Psi Bank 2:")
# print(psi2[0][-1])
# print(psi2[1][-1])
# print("Final Psi Bank 3:")
# print(psi3[0][-1])
# print(psi3[1][-1])
# print("Final Psi Bank 4:")
# print(psi4[0][-1])
# print(psi4[1][-1])
# print("Final Psi Bank 5:")
# print(psi5[0][-1])
# print(psi5[1][-1])
# print("Final Psi Bank 6:")
# print(psi6[0][-1])
# print(psi6[1][-1])

tnew=np.reshape(t, ((k+1),1))
assetDummy=0
trows= len(t)
while (assetDummy < assets):
	print("--Final q Asset " + repr(assetDummy + 1) + ":")
	print(q.item((assetDummy,k)))
	print "Price Impact--"
	print (plot_ft[-1]-q.item((assetDummy,k)))
	assetDummy=assetDummy+1

plt.figure(1)
plt.plot(tnew, allTheta[0], color='orange') 
plt.plot(tnew, allTheta[1], color='blue') 
plt.plot(tnew, allTheta[2], color='yellow') 
plt.plot(tnew, allTheta[3], color='brown') 
plt.plot(tnew, allTheta[4], color='black') 
plt.plot(tnew, allTheta[5], color='cyan') 
plt.plot(t,-ratio1, color='red')  
plt.plot(t,-ratio2, color='green')  	
plt.plot(t,-ratio3, color='maroon')  
plt.plot(t,-ratio4, color='grey')  
plt.plot(t,-ratio5, color='lime')  
plt.plot(t,-ratio6, color='magenta')  	
counter=1
# for row in q:
# 	if counter==1:
# 		plt.plot(t, row, color='plum')
# 	if counter==2:
# 		plt.plot(t, row, color='purple')
# 	else:
# 		plt.plot(t, row, color='lavender')
# 	counter=counter+1
plt.title("Greece Stress Test Results")
plt.xlabel("Time (t)")
plt.ylabel("Normalized Simulation Value")
plum1_patch = mpatches.Patch(color='plum', label='Asset 1 Price')
plum2_patch = mpatches.Patch(color='purple', label='Asset 2 Price')
# lavender_patch = mpatches.Patch(color='lavender', label='Asset 3 Price')
orange_patch = mpatches.Patch(color='orange', label='Capital Ratio Bank 1')
blue_patch = mpatches.Patch(color='blue', label='Capital Ratio Bank 2')
yellow_patch = mpatches.Patch(color='yellow', label='Capital Ratio Bank 3')
brown_patch = mpatches.Patch(color='brown', label='Capital Ratio Bank 4')
black_patch = mpatches.Patch(color='black', label='Capital Ratio Bank 5')
cyan_patch = mpatches.Patch(color='cyan', label='Capital Ratio Bank 6')
red_patch = mpatches.Patch(color='red', label='Selling Ratio Bank 1')
green_patch = mpatches.Patch(color='green', label='Selling Ratio Bank 2')
maroon_patch = mpatches.Patch(color='maroon', label='Selling Ratio Bank 3')
grey_patch = mpatches.Patch(color='grey', label='Selling Ratio Bank 4')
lime_patch = mpatches.Patch(color='lime', label='Selling Ratio Bank 5')
magenta_patch = mpatches.Patch(color='magenta', label='Selling Ratio Bank 6')
plt.legend(handles=[orange_patch, blue_patch, red_patch, green_patch, plum1_patch, plum2_patch])
plt.show()
plt.figure(3)  		
counter=1
for row in q:
	if counter==1:
		plt.plot(t, row, color='plum')
	if counter==2:
		plt.plot(t, row, color='purple')	
	if counter==3:
		plt.plot(t, row, color='red')
	if counter==4:
		plt.plot(t, row, color='blue')
	if counter==5:
		plt.plot(t, row, color='black')	
	if counter==6:
		plt.plot(t, row, color='orange')
	counter=counter+1
plt.plot(plot_t, plot_ft, color='green')
plt.title("Price Impact on Asset Value with D = 81.04 ")
plt.xlabel("Time (t)")
plt.ylabel("Normalized Simulation Value")
plum_patch = mpatches.Patch(color='plum', label='Asset 1 Price')
red_patch = mpatches.Patch(color='red', label='Asset 3 Price')
green_patch = mpatches.Patch(color='green', label='Asset Decay Without Price Impact')
purple_patch = mpatches.Patch(color='purple', label='Asset 2 Price')
blue_patch = mpatches.Patch(color='blue', label='Asset 4 Price')
black_patch = mpatches.Patch(color='black', label='Asset 5 Price')
orange_patch = mpatches.Patch(color='orange', label='Asset 6 Price')
plt.legend(handles=[plum_patch, purple_patch, red_patch, blue_patch, black_patch,orange_patch, green_patch])
plt.show()
