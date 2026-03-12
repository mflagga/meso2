import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

dane = np.loadtxt('psi.dat')
misc = np.loadtxt('misc.dat')
V = np.loadtxt('V.dat')

xmin = misc[0]
xmax = misc[1]
Vmax = misc[2]
bar = misc[3]

tu = np.unique(dane[:,0])

cwd = Path.cwd()
frame_dir = cwd / "frames"
frame_dir.mkdir(parents=True, exist_ok=True)

plt.figure(figsize=(8,6))

n=0
psimax = max(dane[:,4])
psimin = min(dane[:,2])

if (bar):
    ymax = 1.2*max(psimax,Vmax)
else:
    ymax = 1.2*psimax

for t in tu:
    mask = dane[:,0]==t
    plt.plot(V[:,0],V[:,1],color='black',label='V')
    plt.plot(dane[mask,1],dane[mask,2],label='Re',color='blue',linewidth=0.5)
    plt.plot(dane[mask,1],dane[mask,3],label='Im',color='red',linewidth=0.5)
    plt.plot(dane[mask,1],dane[mask,4],label='Norm',color='green')
    #plt.tight_layout()
    plt.title(f't={t:.3f}')
    plt.axvline(x=xmin,ymin=psimin,ymax=psimax,color='black')
    plt.axvline(x=xmax,ymin=psimin,ymax=psimax,color='black')
    plt.ylim(psimin,psimax)
    plt.xlim(1.5*xmin,1.5*xmax)
    plt.xlabel(rf'$x$')
    plt.legend()
    plt.ylabel(rf'$\Psi$')
    #plt.grid(ls=":")
    filename = frame_dir / f"frame_{n:04d}.png"
    n+=1
    plt.savefig(filename)
    plt.clf()
plt.close()