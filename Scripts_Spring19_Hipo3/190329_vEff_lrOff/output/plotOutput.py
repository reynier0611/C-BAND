from __future__ import division
import numpy as np
import sys
import matplotlib.pyplot as plt

fadc_short = []
tdc_short = []
fadc_long = []
tdc_long = []
with open("TEST_effective_velocity.txt") as g:
	for line in g:
		li = line.strip().split("\t")

		# sector | layer | component | veff_tdc | veff_fadc | err | err
		#   0    |   1   |     2     |     3    |     4     |  5  |  6

		if( float(li[3]) == 0 or float(li[4] == 0) ): continue
		if( int(li[0]) == 3 or int(li[0]) == 4 ):
			tdc_short.append( float(li[3] ))
			fadc_short.append( float(li[4] ))
		else:
			tdc_long.append( float(li[3] ))
			fadc_long.append( float(li[4] ))

plt.scatter(fadc_short,tdc_short,color='red',label='Short Bars')
plt.scatter(fadc_long,tdc_long,color='blue',label='Long Bars')
plt.legend(numpoints=1,loc=4)
plt.xlabel('Speed of Light from FADC [cm/ns]',fontsize=14)
plt.ylabel('Speed of Light from TDC [cm/ns]',fontsize=14)
plt.xlim([11,15.5])
plt.ylim([11,15.5])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig("veff.pdf")
plt.show()
