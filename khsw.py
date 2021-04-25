import numpy as np
import matplotlib.pyplot as plt
import sys


arg = sys.argv
M = int(arg[1])

w1 = np.loadtxt('spec{}sub.txt'.format(M))



for i in range(1,M):
    plt.plot(w1[:,0],w1[:,i], color='orange', lw=0.1)

plt.plot(w1[:,0],w1[:,M],label='$\it{SWenergy}$', color='orange', lw=0.1)

plt.legend(fontsize=20)



#plt.xlabel('kx')
plt.ylabel('$\it{ω}$ ',fontsize=27)
plt.xlim([0,8/3])
plt.ylim([0,4])
plt.xticks([0, 1.33, 1.66,8/3],["$\it{Γ}$","$\it{K}$","$\it{M}$","$\it{Γ}$"],fontsize=27)
plt.yticks([0, 1, 2, 3, 4],fontsize=27)


plt.vlines(0, 0, 4, linestyle='dashed', linewidth=0.3,) 
plt.vlines(1.3333, 0, 4, linestyle='dashed', linewidth=0.3) 
plt.vlines(1.6666, 0, 4, linestyle='dashed', linewidth=0.3)
plt.vlines(1, 0, 4, linestyle='dashed',color='red', linewidth=0.8)
plt.vlines(2.1666, 0, 4, linestyle='dashed', color='green', linewidth=0.8)
 




plt.savefig('khsw.pdf',transparent=True,bbox_inches='tight')