from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

rey_tdc = []
efrain_tdc = []
ids = []
with open("../TDC_layer_offsets.txt","rb") as rey, open("TDC_layer_offsets.txt","rb") as efrain:
	for line in rey:
		parse = line.strip().split("\t")
		ids.append( int(parse[1])*100 + int(parse[0])*10 + int(parse[2]) )
		rey_tdc.append( parse[3] )
	for line in efrain:
		parse = line.strip().split("\t")
		efrain_tdc.append( parse[3] )

ids = np.asarray(ids,dtype=int)
rey_tdc = np.asarray(rey_tdc,dtype=float)
efrain_tdc = np.asarray(efrain_tdc,dtype=float)

plt.scatter(ids,rey_tdc-efrain_tdc,marker='o')
plt.xlabel("Bar ID [layer,sector,comp]")
plt.ylabel("(Rey's - Efrain's) TDC Paddle Offsets")
plt.show()
