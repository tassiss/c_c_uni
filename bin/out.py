import glob
import numpy as np
import matplotlib.pyplot as plt
f= glob.glob('../output/*.dat')
f.sort()
l=len(f)
for i in range (0,l):
	data=np.loadtxt(f[i])
	t,x,temp=zip(*data)
	#print data
	#print f[i]
	plt.title('t: '+str(t[0]))
	plt.xlabel('Posicao do ponto')
	plt.ylabel('temperatura')
	plt.plot(x,temp)
	plt.savefig('../output/'+str(t[0])+'.png')
	#plt.show()	
	i=i+1
	
