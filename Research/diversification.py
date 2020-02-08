#algorithm to provide initial asset ownership values for our model based on varying degrees of portfolio diversity.

import numpy as np
Banks = 6
Assets = 6



bank1 =np.array([0,0,0, 0, 0, 0])
bank2 =np.array([0,0,0, 0, 0, 0])
bank3 =np.array([0,0,0, 0, 0, 0])
bank4 =np.array([0,0,0, 0, 0, 0])
bank5 =np.array([0,0,0, 0, 0, 0])
bank6 =np.array([0,0,0, 0, 0, 0])
allBanks = np.array([bank1, bank2, bank3, bank4, bank5, bank6])
allAssets = np.array([100, 100, 100, 100, 100, 100])


degree = 0
while (degree < 37):
	randomAsset = np.random.choice(Assets,1)
	print allAssets[randomAsset-1]
	randomOwnershipProportion = np.random.choice(allAssets[randomAsset-1] ,1)
	print randomOwnershipProportion
	randomBank = np.random.choice(Banks,1)
	if allBanks[randomBank-1,randomAsset-1] == 0:
		allBanks[randomBank-1,randomAsset-1]= randomOwnershipProportion
		allAssets[randomAsset-1] = allAssets[randomAsset-1] - randomOwnershipProportion
		print degree
		print allBanks
		degree = degree +1

