import numpy as np, nibabel as nib, scipy.io as sio, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
D='lab3-for_students/'; O='figs/'
a=nib.load(D+'collab_dwi.nii.gz').get_fdata()[:,:,0,:]
m=np.asarray(nib.load(D+'collab_mask.nii.gz').dataobj).reshape(64,64)>0
b=np.loadtxt(D+'collab_dwi.bval'); g=np.loadtxt(D+'collab_dwi.bvec').T
mo=sio.loadmat('octave_out/map_orig.mat'); mf=sio.loadmat('octave_out/map_fixed.mat')
# voxels in qMRLab/MATLAB 1-based (i,j) -> python (i-1,j-1); image shown as rows=j (y), cols=i (x)
VOX={'eau (12,20)':(12,20),'tube horizontal (16,32)':(16,32),'tube vertical (32,48)':(32,48),'croisement (32,32)':(32,32)}
def show(ax,img,cmap,vmin,vmax,title,marks=True):
    im=ax.imshow(np.where(m,img,np.nan).T,origin='lower',cmap=cmap,vmin=vmin,vmax=vmax)
    ax.set_title(title); ax.axis('off')
    if marks:
        for n,(i,j) in VOX.items(): ax.plot(i-1,j-1,'s',ms=6,mfc='none',mec='k',mew=1.2)
    plt.colorbar(im,ax=ax,fraction=0.046,pad=0.03)
# ---------- Fig 3A : reproduction
fig,ax=plt.subplots(1,4,figsize=(15,3.6))
show(ax[0],mo['FA'],'viridis',0,1.2,'FA (original)')
for k in range(3): show(ax[k+1],mo['L'][...,k],'RdBu_r',-3,3,f'L{k+1} (original), µm²/ms')
fig.suptitle('Fig. 3A – Fit qMRLab (dti, linéaire) des données telles que reçues. Carrés : voxels inspectés',fontsize=10)
plt.tight_layout(); plt.savefig(O+'fig3A_maps_original.png',dpi=150); plt.close()
rows=[]
for n,(i,j) in VOX.items():
    rows.append((n,*[mo['L'][i-1,j-1,k] for k in range(3)],mo['FA'][i-1,j-1],*[mf['L'][i-1,j-1,k] for k in range(3)],mf['FA'][i-1,j-1]))
np.save(O+'vox_table.npy',np.array(rows,dtype=object))
for r in rows: print(r[0],' orig L=%.2f %.2f %.2f FA=%.2f | fixed L=%.2f %.2f %.2f FA=%.2f'%tuple(r[1:]))
# ---------- Fig 3B-1 : raw volumes with labels
fig,ax=plt.subplots(3,7,figsize=(15,7))
lab={(1,0,0):'x',(0,1,0):'y',(0,0,1):'z',(1,1,0):'xy',(1,0,1):'xz',(0,1,1):'yz',(0,0,0):'–'}
for k in range(19):
    A=ax.flat[k]; A.imshow(a[...,k].T,origin='lower',cmap='gray',vmin=0,vmax=1000)
    A.set_title(f'vol {k+1}: b={b[k]:.0f}, dir {lab[tuple((g[k]>0.1).astype(int))]}',fontsize=8,color='C3' if k in (0,18) else 'k'); A.axis('off')
for A in ax.flat[19:]: A.axis('off')
fig.suptitle('Fig. 3B-1 – Les 19 volumes bruts avec l\'étiquette du protocole (même fenêtrage 0–1000). Le volume étiqueté b=0 (vol 1) est plus sombre que le volume étiqueté b=1500 (vol 19), qui est le plus brillant de tous.',fontsize=9)
plt.tight_layout(); plt.savefig(O+'fig3B1_raw_volumes.png',dpi=130); plt.close()
# ---------- Fig 3B-2 : water ROI signal vs labelled b
water=m.copy(); water[22:42,:]=False; water[:,22:42]=False
S=a[water].mean(0); S0=S[0]
fig,ax=plt.subplots(1,2,figsize=(11,3.8))
dirs=['x','y','z','xy','xz','yz']
for d in range(6):
    idx=[0]+list(range(1+3*d,4+3*d))
    ax[0].plot(b[idx],S[idx]/S0,'o-',lw=1.5,ms=5,label=dirs[d])
bb=np.linspace(0,1500,50); ax[0].plot(bb,np.exp(-bb/1000*2),'k--',lw=1.5,label='attendu (D = 2 µm²/ms)')
ax[0].set_xlabel('b étiqueté (s/mm²)'); ax[0].set_ylabel('S / S(vol 1, « b=0 »)'); ax[0].set_title('Avec l\'ordre reçu'); ax[0].legend(fontsize=7,ncol=2); ax[0].grid(alpha=.3)
Sr=S[::-1]
for d in range(6):
    idx=[0]+list(range(1+3*d,4+3*d))
    ax[1].plot(b[idx],Sr[idx]/Sr[0],'o-',lw=1.5,ms=5,label=dirs[d])
ax[1].plot(bb,np.exp(-bb/1000*2),'k--',lw=1.5,label='attendu (D = 2 µm²/ms)'); ax[1].set_yscale('log')
ax[1].set_xlabel('b (s/mm²)'); ax[1].set_ylabel('S / S0 (échelle log)'); ax[1].set_title('Après inversion de l\'ordre des volumes'); ax[1].legend(fontsize=7,ncol=2); ax[1].grid(alpha=.3)
fig.suptitle(f'Fig. 3B-2 – Signal moyen dans l\'eau ({water.sum()} voxels, hors tubes) en fonction du b. À gauche le signal AUGMENTE avec b (S/S0 > 1) : physiquement impossible.',fontsize=9)
plt.tight_layout(); plt.savefig(O+'fig3B2_water_signal_vs_b.png',dpi=150); plt.close()
print('water ratio orig', (S/S0).round(2))

# ===================== partie 2 =====================
import numpy as np, nibabel as nib, scipy.io as sio, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
D='lab3-for_students/'; O='figs/'
a=nib.load(D+'collab_dwi.nii.gz').get_fdata()[:,:,0,:]
m=np.asarray(nib.load(D+'collab_mask.nii.gz').dataobj).reshape(64,64)>0
b=np.loadtxt(D+'collab_dwi.bval')/1000; g=np.loadtxt(D+'collab_dwi.bvec').T
mo=sio.loadmat('octave_out/map_orig.mat'); mf=sio.loadmat('octave_out/map_fixed.mat'); sw=sio.loadmat('octave_out/sim_water.mat')
cols={0.5:'#2a6fdb',1.0:'#e08a00',1.5:'#c8326b',0:'#555555'}
def qplot(ax,S,Dv,title):
    Dt=np.asarray(Dv,float).reshape(3,3,order='F'); w,V=np.linalg.eigh(Dt); f=V[:,np.argmax(w)]
    S0=S[b==0].mean(); x=g@f
    for bv in [0,0.5,1.0,1.5]:
        k=np.where(np.isclose(b,bv))[0]; o=np.argsort(x[k]); k=k[o]
        ax.plot(x[k],S[k],'x',color=cols[bv],ms=7,mew=1.6,label=f'b={bv*1000:.0f}')
        model=S0*np.exp(-bv*np.einsum('ij,jk,ik->i',g[k],Dt,g[k]))
        ax.plot(x[k],model,'-',color=cols[bv],lw=1.5)
    ax.set_xlabel('G∥/|G|  (projection sur la direction principale)'); ax.set_ylabel('Signal'); ax.set_title(title,fontsize=9); ax.grid(alpha=.3); ax.legend(fontsize=7)
i,j=12,20
fig,ax=plt.subplots(1,3,figsize=(15,4))
L=mo['L'][i-1,j-1]; qplot(ax[0],a[i-1,j-1],mo['Dm'][i-1,j-1],f'Voxel eau ({i},{j}), données reçues\nL = {L[0]:.2f}, {L[1]:.2f}, {L[2]:.2f} µm²/ms ; FA = {mo["FA"][i-1,j-1]:.2f}')
Ls=sw['Lsim'].ravel(); qplot(ax[1],sw['Ssim'].ravel(),sw['Dsim'].ravel(),f'Simulation qMRLab : eau L = 2, 2, 2 ; SNR = 50\nfit : L = {Ls[0]:.2f}, {Ls[1]:.2f}, {Ls[2]:.2f} ; FA = {Ls[3]:.2f}')
L=mf['L'][i-1,j-1]; qplot(ax[2],a[i-1,j-1,::-1],mf['Dm'][i-1,j-1],f'Même voxel, volumes remis dans l\'ordre\nL = {L[0]:.2f}, {L[1]:.2f}, {L[2]:.2f} µm²/ms ; FA = {mf["FA"][i-1,j-1]:.2f}')
fig.suptitle('Fig. 3B-3 – Croix = données, traits = modèle. Données reçues : le signal « b=1500 » est le plus haut et le « b=0 » est au niveau des b=1500 ; la simulation de l\'eau attendue montre l\'inverse.',fontsize=9)
plt.tight_layout(); plt.savefig(O+'fig3B3_voxelfit_vs_simulation.png',dpi=150); plt.close()
# ---- 3E maps after repair
VOX={(12,20),(16,32),(32,48),(32,32)}
def show(ax,img,cmap,vmin,vmax,title):
    im=ax.imshow(np.where(m,img,np.nan).T,origin='lower',cmap=cmap,vmin=vmin,vmax=vmax); ax.set_title(title); ax.axis('off')
    for (i,j) in VOX: ax.plot(i-1,j-1,'s',ms=6,mfc='none',mec='w' if cmap=='viridis' else 'k',mew=1.2)
    plt.colorbar(im,ax=ax,fraction=0.046,pad=0.03)
fig,ax=plt.subplots(1,4,figsize=(15,3.6))
show(ax[0],mf['FA'],'viridis',0,1,'FA (réparé)')
for k in range(3): show(ax[k+1],mf['L'][...,k],'magma',0,2.5,f'L{k+1} (réparé), µm²/ms')
fig.suptitle('Fig. 3E-1 – Fit qMRLab après réparation (volumes remis dans l\'ordre). Eau : L ≈ 2 µm²/ms partout, FA ≈ 0.1 ; tubes : cigare ; croisement : FA intermédiaire.',fontsize=9)
plt.tight_layout(); plt.savefig(O+'fig3E1_maps_fixed.png',dpi=150); plt.close()
# ---- 3E principal direction
fig,ax=plt.subplots(1,2,figsize=(10,4.2))
for A,mp,t in [(ax[0],mo,'Données reçues'),(ax[1],mf,'Après réparation')]:
    Dt=mp['Dm'].reshape(64,64,3,3,order='F'); w,V=np.linalg.eigh(Dt); v1=V[...,:,-1]
    rgb=np.abs(v1)*np.clip(mp['FA'],0,1)[...,None]; rgb[~m]=1
    A.imshow(np.transpose(rgb,(1,0,2)),origin='lower'); A.axis('off')
    A.set_title(t+'\nCouleur = |V1| (rouge x, vert y, bleu z) × FA',fontsize=9)
fig.suptitle('Fig. 3E-2 – Direction principale de diffusion. Après réparation : tube horizontal en rouge (x), tube vertical en vert (y), eau quasi noire (FA≈0).',fontsize=9)
plt.tight_layout(); plt.savefig(O+'fig3E2_V1_direction.png',dpi=150); plt.close()
# ROI summary
water=m.copy(); water[22:42,:]=False; water[:,22:42]=False
th=np.zeros_like(m); th[6:26,28:35]=True; th[38:58,28:35]=True; th&=m
tv=np.zeros_like(m); tv[28:35,6:26]=True; tv[28:35,38:58]=True; tv&=m
cr=np.zeros_like(m); cr[29:35,29:35]=True
for n,R in [('orig',mo),('fixed',mf)]:
  for rn,M in [('eau',water),('tube x',th),('tube y',tv),('croisement',cr)]:
    L=R['L'][M]; print(n,rn,M.sum(),'L1 %.2f±%.2f L2 %.2f±%.2f L3 %.2f±%.2f FA %.2f±%.2f'%(L[:,0].mean(),L[:,0].std(),L[:,1].mean(),L[:,1].std(),L[:,2].mean(),L[:,2].std(),R['FA'][M].mean(),R['FA'][M].std()))
