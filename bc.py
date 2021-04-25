import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as matcol
from matplotlib.cm import ScalarMappable, get_cmap
import sys

texparams = {'ps.useafm': True, 'pdf.use14corefonts': True, 'pdf.fonttype': 42,
             'text.usetex': True, 'text.latex.preamble': r'\usepackage{amsmath} \usepackage{txfonts} \usepackage{bm}'}
plt.rcParams.update(texparams)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

norm = matcol.Normalize(vmin=-1, vmax=1)
cmap = plt.get_cmap("coolwarm")

def cplot(ax,x,y,z):
    c = cmap(norm(z)) 
    for j in range(len(x)-1):
        ax.plot(x[j:j+2], y[j:j+2], color=c[j+1], lw=1)

f = np.loadtxt('spec4sub.txt')
g = np.loadtxt('ch4sub.txt')


fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111)

arg = sys.argv
M = int(arg[1])

for i in range(1,M+1,1)
  G1 = cplot(ax,f[:,0],f[:,M],g[:,M])

#ax.legend(handles=G1,labels=[r"$LSWT$"], fontsize=20)

yposi = np.array([0, 1, 2, 3, 4])
xposi = np.array([0, 1.33, 1.66, 8./3.])
xposi_l = [r'$\Gamma$',r'$K$',r'$M$',r'$\Gamma$']

ax.set_xticks(xposi)
ax.set_xticklabels(xposi_l, fontsize=20)
ax.set_yticks(yposi)
ax.set_yticklabels(yposi,fontsize=20)
#ax.set_ylabel(r'$\omega$',fontsize=24)
ax.set_xlim(xposi[0],xposi[3])
ax.set_ylim(yposi[0],yposi[yposi.size-1])

for x in xposi:
    ax.axvline(x, color='black', ls='dashed', lw=0.5)

mappable = ScalarMappable(cmap=cmap, norm=norm)
mappable._A = []

#plt.xticks(color="None")
plt.yticks(color="None")

#ts = [-2, -1, 0, 1, 2]
ts = [-1,-0.5,0,0.5,1]
# c = fig.colorbar(mappable,ticks=ts)
# c.ax.set_yticklabels(ts, fontsize=20)







plt.savefig('bc.pdf',transparent=True, bbox_inches='tight')
