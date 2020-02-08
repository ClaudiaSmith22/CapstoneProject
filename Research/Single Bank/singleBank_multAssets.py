
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from main import *

#Input Initial Values
b=.04	#constant in f(gamma(t))
s_initial=np.array([[20.0],[11.0]])  #initial stocks
s=np.copy(s_initial) #array of current stock values (used in algorithm)
dt=0.0001  
start=0  #Holds a value of 1 if bank is selling
t=np.array([0.0])
alpha=10.0
theta_min= 0.08
event=np.array([-5])
Gamma= np.array([[0.0], [0.0]]) #Inital number of stocks sold at t. Always 0.
psi= np.array([[0.0], [0.0]]) #Initial cash made from liquidation. Always 0.
q = np.array([[1.0], [1.0]])  #Price of stock
ratio=np.array([0]) #Start out selling 0% of assets

k=0
assets=2 #Insert number of assets

#Itialize F values
f_t= np.exp(-t[k])
plot_ft = np.array([np.exp(-t[k])])
fprime_t= -f_t
f_g=fGamma(-b, Gamma, assets, k)
fprime_g=fGamma_prime(b, f_g, assets)

#Initialize theta
theta_dot=updateTheta2(event, psi, s_initial, ratio.item((k)), q, alpha, k)
theta=np.array(theta_dot)

while (t[k]<1): 
	print t[k]

	#update F values
	f_t= np.exp(-t[k])
	plot_ft = np.concatenate((plot_ft, [np.exp(-t[k])]))
	fprime_t= -f_t
	f_g=fGamma(b, Gamma, assets, k)
	fprime_g=fGamma_prime(b, f_g, assets)

	#calculate gamma, psi, t, z, q with previous k values
	n=0 
	for row in Gamma:
		Gamma_new=updateGamma(start, s[n][-1], ratio.item((k)))
		if n==0:
			Gamma_dot=np.array([Gamma.item((n,k))+Gamma_new*dt])
			Psi_dot=np.array([psi.item((n,k))+Gamma_new*q.item((n,k))*dt])
		else:
			Gamma_dot=np.concatenate((Gamma_dot, [Gamma.item((n,k))+Gamma_new*dt]))
			Psi_dot=np.concatenate((Psi_dot, [psi.item((n,k))+Gamma_new*q.item((n,k))*dt]))
		n=n+1

	#update q values
	q_dot=f_t*f_g
	#Reshape arrays and append new values
	Gamma_dot=(np.reshape(Gamma_dot, (n,1))) 
	Psi_dot=(np.reshape(Psi_dot, (n,1)))
	q = np.concatenate((q, np.reshape(q_dot, (n,1))), axis=1)
	Gamma= np.concatenate((Gamma, Gamma_dot), axis=1)
	psi = np.concatenate((psi, Psi_dot), axis=1)
	t= np.concatenate((t, [t[k]+dt]))

	#Calculate Theta last
	theta_dot=updateTheta2(event, psi, s_initial, ratio.item((k)), q, alpha, k)
	theta=np.concatenate((theta, theta_dot))

	#Update current stock units
	s_dot= np.array(np.reshape(s[:,0], (n,1))-Gamma_dot)
	s= np.concatenate((s, np.reshape(s_dot, (n,1))), axis=1)

	#Update k
	k=k+1
	
	#Check if below the threshold
	if (theta[k]-theta_min < 1e-20):
		#Set the start value to 1 to indicate selling, find ratio
		start=1
		ratio=np.concatenate((ratio, [findRatio(event, psi, q[:,k], theta_min, alpha, s_initial, k)]))
		counter=2
		#Delete the theta, gamma, psi, and q values for the last event t
		theta=np.delete(theta, k, 0)
		Gamma=np.delete(Gamma, k, 1)
		psi=np.delete(psi, k, 1)
		q=np.delete(q, k, 1)
		plot_ft=np.delete(plot_ft, k, 0)
		limit=t[k]
		t=np.delete(t, k, None)
		
		#Set the dt value to 100th of the original
		dt_new=dt/100

		#Update the itertion count to the previous
		k=k-1
		c=0

		while (limit-t[k]>dt/100):  #iterate until "whole" interval 
			#Update F values
			f_t= np.exp(-t[k])
			plot_ft = np.concatenate((plot_ft, [np.exp(-t[k])]))
			fprime_t= -f_t
			f_g=fGamma(-b, Gamma, assets, k)
			fprime_g=fGamma_prime(b, f_g, assets)

			#Update z, gamma, q, psi, t and theta
			n=0 
			for row in Gamma:
				Gamma_new=updateGamma(start, s[n][-1], ratio.item((k)))
				if n==0:
					Gamma_dot=np.array([Gamma.item((n,k))+Gamma_new*dt_new])
					Psi_dot=np.array([psi.item((n,k))+Gamma_new*q.item((n,k))*dt_new])
				else:
					Gamma_dot=np.concatenate((Gamma_dot, [Gamma.item((n,k))+Gamma_new*dt_new]))
					Psi_dot=np.concatenate((Psi_dot, [psi.item((n,k))+Gamma_new*q.item((n,k))*dt_new]))
				n=n+1

			#Update q values
			q_dot=f_t*f_g

			#Reshape arrays and append new values
			Gamma_dot=(np.reshape(Gamma_dot, (n,1)))
			Psi_dot=(np.reshape(Psi_dot, (n,1)))
			q = np.concatenate((q, np.reshape(q_dot, (n,1))), axis=1)
			Gamma= np.concatenate((Gamma, Gamma_dot), axis=1)
			psi = np.concatenate((psi, Psi_dot), axis=1)
			t=np.concatenate((t, [t[k]+dt_new]))
			theta_dot=updateTheta2(event, psi, s_initial, ratio.item(k), q, alpha, k)
			theta=np.concatenate((theta, theta_dot))

			#Update current stock units
			s_dot= np.array(np.reshape(s[:,0], (n,1))-Gamma_dot)
			s= np.concatenate((s, np.reshape(s_dot, (n,1))), axis=1)
			if c>0:
				ratio=np.concatenate((ratio, [ratio.item((k))]))		
			c=c+1
			k=k+1
	else:
		ratio=np.concatenate((ratio, [ratio.item((k-1))]))



#Print final values and graphs
print("Final Theta:")
print (theta[k])
print("Final t:")
print(t[k])
print("Final S asset 1:")
print (s[0][-1])
print("Final S asset 2:")
print(s[1][-1])
tnew=np.reshape(t, ((k+1),1))
assetDummy=0
trows= len(t)
while (assetDummy < assets):
	print("Final Gamma Asset " + repr(assetDummy + 1) + ":")
	print(Gamma.item((assetDummy,k)))
	print("Final q Asset " + repr(assetDummy + 1) + ":")
	print(q.item((assetDummy,k)))
	print("Final Psi Asset " + repr(assetDummy + 1) + ":")
	print(psi.item((assetDummy,k)))
	assetDummy=assetDummy+1
plt.figure(1)
plt.plot(tnew, theta) 
plt.plot(t,-ratio)  		
for row in q:
	plt.plot(t, row)
plt.title("Single-Bank, Multi-Asset Method 1 Results")
plt.xlabel("Time (t)")
plt.ylabel("Scaled Simulation Values")
red_patch = mpatches.Patch(color='red', label='Asset 1 Price')
green_patch = mpatches.Patch(color='green', label='Asset 2 Price')
blue_patch = mpatches.Patch(color='blue', label='Capital Ratio')
orange_patch = mpatches.Patch(color='orange', label='Selling Ratio')
plt.legend(handles=[red_patch, green_patch, blue_patch, orange_patch])
plt.show()
plt.figure(2)
plt.plot(tnew, theta) 
plt.plot(t,-ratio)  		
plt.title("Single-Bank, Multi-Asset Method 1 Selling and Capital Ratio Results")
plt.xlabel("Time (t)")
plt.ylabel("Scaled Simulation Value")
plt.legend(handles=[blue_patch, orange_patch])
plt.show()
plt.figure(3)  		
for row in q:
	plt.plot(t, row)
plt.plot(t, plot_ft)
plt.title("Price Impact on Asset Value")
plt.xlabel("Time (t)")
plt.ylabel("Scaled Simulation Value")
red_patch = mpatches.Patch(color='red', label='Asset 1 Price')
green_patch = mpatches.Patch(color='green', label='Asset Decay Without Price Impact')
blue_patch = mpatches.Patch(color='blue', label='Asset 2 Price')
plt.legend(handles=[red_patch, blue_patch, green_patch])
plt.show()