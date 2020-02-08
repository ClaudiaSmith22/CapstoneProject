import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA

#Functions to be called from the single bank model scripts

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
			sqArray=np.array([s[i]*q[i]])
		else:
			sqArray=np.concatenate((sqArray,[s[i]*q[i]]))
		i=i+1
	psiSum=psiSum.item(((i-1),0))
	numerator=event+psiSum
	sqSum=np.cumsum(sqArray, axis=0)
	sqSum=sqSum[i-1]
	denominator=sqSum*(theta_min*alpha-1)
	ratio=1-(numerator/denominator)
	ratio=ratio[0]
	return ratio

def fGamma(b, Gamma, assets, k):
	i=0
	while (i<assets):
		if (i==0):
			f_g=np.array([np.exp(-b*Gamma.item(i,k))])
		else:
			f_g=np.concatenate((f_g, [np.exp(-b*Gamma.item(i,k))]))
		i=i+1
	return f_g

def fGamma_prime(b, f_g, assets):
	j=0
	while (j<assets):
		if (j==0):
			fprime_g=np.array([-b*f_g[j]])
		else:
			fprime_g=np.concatenate((fprime_g,[-b*f_g[j]]))
		j=j+1
	return fprime_g

def updateTheta2(event, psi, s, ratio, q, alpha, k):
	n=0
	for row in psi:
		sum1=psi.item(n,k)+s[n]*(1-ratio)*q.item(n,k)
		sum2=alpha*(s[n]*(1-ratio)*q.item(n,k))
		if sum1<0:
			sum1=[0]
		if n==0:
			theta1=np.array([sum1])
			theta2=np.array([sum2])
		else:
			theta1=np.concatenate((theta1, [sum1]))
			theta2=np.concatenate((theta2, [sum2]))
		n=n+1
	numSum= np.cumsum(theta1, axis=0)
	denSum= np.cumsum(theta2, axis=0)
	numerator= event+numSum.item(((n-1),0))
	if numerator <0 :
		numerator=0
	denominator=denSum.item(((n-1),0))
	theta_dot=numerator/denominator
	theta_dot=np.reshape(theta_dot, ((n-1),1))
	return theta_dot


def findRatioDeriv(event, psi, q, s, k, q_dot, ratio, dt):
	psiSum=np.cumsum(psi, axis=0)
	psi_dotSum=np.cumsum(psi_Dot, axis=0)
	i=0
	for row in psi:
		if i==0:
			sqArray=np.array([s.item((i,k))*q[i]])
		else:
			sqArray=np.concatenate((sqArray,[s.item((i,k))*q[i]]))
		i=i+1

	i=0
	for row in psi:
		if i==0:
			sDqArray=np.array([s.item((i,k))*q_dot[i]])
		else:
			sDqArray=np.concatenate((sDqArray,[s.item((i,k))*q_dot[i]]))
		i=i+1

	psiSum=psiSum.item(((i-1),0))
	psi_dotSum=psi_dotSum.item((i-1),0)
	sqSum=np.cumsum(sqArray, axis=0)
	sqSum=sqSum[i-1]
	sDqSum=np.cumsum(sDqArray, axis=0)
	sDqSum=sDqSum[i-1]

	one=((1-ratio.item(k-1)) * psi_dotSum) / (event+psiSum)
	two = sDqSum/sqSum
	three = ((ratio.item(k-1))*sDqSum) / sqSum
	ratio= ratio.item(k-1) + one-two+three

	ratio=ratio[0]*dt
	return ratio






















# def updateZ(event, psi, s, Gamma, q):
# 	z1=(-event-psi)*(s-Gamma)
# 	z2=q*(q*(s-Gamma)+psi+event)
# 	if z2==0:
# 		z=0
# 	else:
# 		z=z1/z2
# 	return z

# def updateGamma2(z, fprime_t, fprime_g, f_t, f_g, start):
# 	Gamma1=z*fprime_t*f_g
# 	Gamma2=1+(z*f_t*fprime_g)
# 	if (start==1):
# 		Gamma_dot=-Gamma1/Gamma2
# 	else:
# 		Gamma_dot= 0
# 	return Gamma_dot

# def updateTheta(event, psi, s, Gamma, q, alpha):
# 	theta1=event+psi+(s-Gamma)*q
# 	if theta1<0:
# 		theta1=0
# 	theta2=alpha*(s-Gamma)*q 
# 	theta=theta1/theta2
# 	return theta