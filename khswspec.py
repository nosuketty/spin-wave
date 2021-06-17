from locale import YESEXPR
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import sys


from matplotlib.colors import LinearSegmentedColormap

arg = sys.argv
M = int(arg[1])

def generate_cmap(colors):
    """自分で定義したカラーマップを返す"""
    values = range(len(colors))

    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append((v / vmax, c))
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)


def generate_cmapv(colors):
    """自分で定義したカラーマップを返す"""
    values = [v[1] for v in colors]
    vmax = np.max(values)
    vmin = np.min(values)
    color_list = []
    for c in colors:
        color_list.append(((c[1]-vmin)/(vmax-vmin), c[0]))
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)


cm = generate_cmapv([('black', 0.), ('navy', 0.1), ('red', 0.2), ('orange', 0.3),
                     ('yellow', 0.6), ('#bbeebb', 0.8), ('#e0ffff', 0.9), ('white', 1)])


def fl(x, x0, gamma):
    return (1./np.pi)*gamma/((x-x0)*(x-x0) + gamma*gamma)


gamma = 0.1


#マイデータ 0=x, 1~8=ene0~7, 9~16=xx, 17~24=yy, 25~32=zz, 33~40=xy, 41~48=xz, 49~56=yz
ndata = np.loadtxt('spec{}sub.txt'.format(M))


X = ndata[:, 0]
Y = np.linspace(-2, 4, 2000)
Z = np.empty((Y.size, X.size))

nw = np.empty((M, Y.size, X.size))
nwa = np.empty((M, X.size))
nweight = np.empty((M, Y.size, X.size))

x, y = np.meshgrid(X, Y)

for i in range(M):
    nw[i,:,:], NULL = np.meshgrid(ndata[:,i+1], Y)

#( ( xx+yy+4*zz+2*xy-4*yz-4*xz )/6 + ( xx+yy-2*xy )/2 )/2
for i in range(0, M):
    nwa[i] = ((ndata[:, M * 1 + 1 + i]+ndata[:, M * 2 + 1 +i]+4*ndata[:, M * 3 + 1 + i]+2*ndata[:, M * 4 + 1 + i]-4*ndata[:, M * 6 + 1 + i]-4*ndata[:, M * 5 + 1 + i])/6.+(ndata[:, M * 1 + 1 + i]+ndata[:, M * 2 + 1 + i]-2*ndata[:, M * 4 + 1 + i])/2.)/2.
    nweight[i], NULL = np.meshgrid(nwa[i], Y)
    Z = Z + nweight[i] * fl(y, nw[i], gamma)

# Z = nweight[0]*fl(y, nw[0], gamma) + nweight[1]*fl(y, nw[1], gamma) + nweight[2]*fl(y, nw[2], gamma) + nweight[3]*fl(y, nw[3], gamma) + nweight[4]*fl(y, nw[4], gamma) + nweight[5]*fl(y, nw[5], gamma) + nweight[6]*fl(y, nw[6], gamma) + nweight[7]*fl(y, nw[7], gamma)

fig = plt.figure()
ax = fig.add_subplot(111)

c = ax.pcolormesh(X, Y, Z, shading='auto', cmap=cm, vmin=0, vmax=3,
                  linewidth=0, rasterized=True)

fig.colorbar(c, ax=ax)

root3 = np.sqrt(3.)

mposi = np.array([0., 1.33, 1.66, 8./3.])
mposi_l = ['$\it{Γ}$','$\it{K}$','$\it{M}$','$\it{Γ}$']
y =np.array([-1,0,1,2,3,4])
y_l = ['-1',  '0', '1', '2', '3', '4']
ax.set_xticks(mposi)
ax.set_xticklabels(mposi_l,fontsize=20)
ax.set_xlim(mposi[0], mposi[3])
ax.set_yticks(y)
ax.set_yticklabels(y_l, fontsize=20)
ax.set_ylim(-2, 4)


#for x in pos:
#    ax.axvline(x, color='w', ls='-', lw=0.5)

for x in mposi:
    ax.axvline(x, color='white', ls='-', lw=0.5)

ax.axhline(y=0, color='white', ls='-', lw=0.5) 

plt.savefig('spec.pdf', transparent=True, bbox_inches='tight')
