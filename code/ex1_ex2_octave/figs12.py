"""Figures exercices 1 et 2 (donnees demo qMRLab dti, fits faits avec le code qMRLab 3.0.0 sous Octave)."""
import numpy as np, nibabel as nib, scipy.io as sio, matplotlib, sys
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
D = 'dti_demo/'; W = 'out/'; O = 'figs/'
Z = 23                     # coupe axiale (0-based) = z 24 en MATLAB
VOX = {'genou': (35, 63, 24), 'cortex': (28, 78, 24), 'ventricule': (34, 54, 24)}   # MATLAB 1-based
NAMES = {'genou': 'genou du corps calleux', 'cortex': 'substance grise corticale', 'ventricule': 'ventricule latéral'}

P = np.loadtxt(W + 'prot.txt')
g = P[:, :3]; G = P[:, 3]; Dl, dl = P[:, 4], P[:, 5]
gamma = 2 * np.pi * 42.576e6
b = (gamma * G * dl) ** 2 * (Dl - dl / 3) / 1e9          # ms/um^2  (= s/mm^2 / 1000)
shell = np.round(b * 1000, -1)
img = nib.load(D + 'DiffusionData.nii.gz')
sims = sio.loadmat(W + 'sims.mat', squeeze_me=True)
COL = {0: '#555555', 300: '#2a6fdb', 1000: '#e08a00', 1990: '#c8326b'}

def qplot(ax, S, Dv, S0, title, true_S=None):
    """Replique de scd_display_qspacedata3D (qMRLab) : x = projection de g sur la direction principale."""
    Dt = np.asarray(Dv, float).reshape(3, 3, order='F'); w, V = np.linalg.eigh(Dt); f = V[:, np.argmax(w)]
    x = g @ f
    for bv in sorted(COL):
        k = np.where(shell == bv)[0]
        if len(k) == 0: continue
        k = k[np.argsort(x[k])]
        if true_S is not None:
            ax.plot(x[k], true_S[k], 'o', mfc='none', color=COL[bv], ms=5, mew=0.8)
        ax.plot(x[k], S[k], 'x', color=COL[bv], ms=5, mew=1.2, label=f'b = {bv:.0f} s/mm²')
        ax.plot(x[k], S0 * np.exp(-b[k] * np.einsum('ij,jk,ik->i', g[k], Dt, g[k])), '-', color=COL[bv], lw=1.4)
    ax.set_xlabel('G∥/|G|'); ax.set_ylabel('Signal'); ax.set_title(title, fontsize=9); ax.grid(alpha=.3)
    ax.legend(fontsize=7, loc='best')

def vol_idx(Gv): return int(np.where(np.isclose(G, Gv, atol=1e-3))[0][0])

def fig1B():
    vols = [(vol_idx(0), 'b = 0'), (vol_idx(0.0566), 'b = 1000'), (vol_idx(0.08), 'b = 1990')]
    sx, sy, sz = 37, 43, Z
    fig, ax = plt.subplots(3, 3, figsize=(10, 10.5))
    lett = iter('abcdefghi')
    for r, (v, lab) in enumerate(vols):
        V = np.asarray(img.dataobj[..., v]).astype(float)
        lo, hi = np.percentile(V[V > 0], [1, 99.5])
        views = [np.rot90(V[:, :, sz]), np.rot90(V[:, sy, :]), np.rot90(V[sx, :, :])]
        for c, (im, nm) in enumerate(zip(views, ['axiale (z = 24)', 'coronale (y = 44)', 'sagittale (x = 38)'])):
            ax[r, c].imshow(im, cmap='gray', vmin=lo, vmax=hi, aspect='equal')
            ax[r, c].set_title(f'Fig. 1B-{next(lett)} – {nm}, {lab} (vol {v + 1})', fontsize=8.5); ax[r, c].axis('off')
    fig.suptitle('Données démo dti (qMRLab) – vues brutes, sans masque, contraste ajusté par volume (1er–99.5e centile)', fontsize=9)
    plt.tight_layout(); plt.savefig(O + 'fig1B_raw.png', dpi=140); plt.close()

def fig1C():
    fa_demo = np.asarray(nib.load(D + 'FitResults/FA.nii.gz').dataobj)[:, :, Z]
    for k, (key, (i, j, z)) in enumerate(VOX.items(), 1):
        fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw={'width_ratios': [1, 1.6]})
        ax[0].imshow(fa_demo.T, origin='lower', cmap='gray', vmin=0, vmax=1)
        ax[0].plot(i - 1, j - 1, '+', color='#e8302a', ms=14, mew=2); ax[0].axis('off')
        ax[0].set_title(f'Voxel [{i}, {j}, {z}] – {NAMES[key]}\n(coupe axiale z = {z}, fond : FA)', fontsize=9)
        L = sims[f'vox{k}_L']
        qplot(ax[1], sims[f'vox{k}_S'], sims[f'vox{k}_D'], sims[f'vox{k}_S0'],
              f'Fit voxel (dti, options par défaut) : FA = {L[3]:.2f} ; L1, L2, L3 = {L[0]:.2f}, {L[1]:.2f}, {L[2]:.2f} µm²/ms')
        fig.suptitle(f'Fig. 1C-{"abc"[k-1]} – Fit d\'un seul voxel : {NAMES[key]}. Croix = données, traits = modèle ajusté.', fontsize=9)
        plt.tight_layout(); plt.savefig(O + f'fig1C_{key}.png', dpi=140); plt.close()

def crop(a, m):
    ii, jj = np.where(m); return a[ii.min() - 2: ii.max() + 3, jj.min() - 2: jj.max() + 3]

def fig1D():
    b0 = np.load(W + 'b0_slice.npy'); m = np.asarray(nib.load(W + 'Mask_slice.nii.gz').dataobj)[:, :, Z] > 0
    fig, ax = plt.subplots(1, 2, figsize=(9, 4.6))
    lo, hi = np.percentile(b0[b0 > 0], [1, 99.5])
    ax[0].imshow(b0.T, origin='lower', cmap='gray', vmin=lo, vmax=hi); ax[0].set_title('b = 0 (moyenne des 9 volumes b = 0), coupe z = 24'); ax[0].axis('off')
    ax[1].imshow(b0.T, origin='lower', cmap='gray', vmin=lo, vmax=hi)
    ov = np.zeros(m.shape + (4,)); ov[m] = [0.9, 0.2, 0.15, 0.35]
    ax[1].imshow(np.transpose(ov, (1, 0, 2)), origin='lower'); ax[1].contour(m.T, levels=[0.5], colors='#e8302a', linewidths=1.2, origin='lower')
    ax[1].set_title(f'Masque superposé ({m.sum()} voxels sur {m.size})'); ax[1].axis('off')
    fig.suptitle('Fig. 1D – Masque du cerveau sur une seule coupe axiale (seuil d\'Otsu + remplissage + plus grande composante + fermeture)', fontsize=9)
    plt.tight_layout(); plt.savefig(O + 'fig1D_mask.png', dpi=140); plt.close()

def fig1E():
    F = sio.loadmat(W + 'fit_slice.mat', squeeze_me=True)
    m = np.asarray(nib.load(W + 'Mask_slice.nii.gz').dataobj)[:, :, Z] > 0
    maps = [('FA', F['FA'], 'viridis', (0, 1), ''), ('L1', F['L'][..., 0], 'magma', (0, 3.5), 'µm²/ms'),
            ('L2', F['L'][..., 1], 'magma', (0, 3.5), 'µm²/ms'), ('L3', F['L'][..., 2], 'magma', (0, 3.5), 'µm²/ms')]
    fig, ax = plt.subplots(1, 4, figsize=(15, 4.6))
    for k, (n, a, cm, lim, u) in enumerate(maps):
        im = ax[k].imshow(crop(np.where(m, a, np.nan), m).T, origin='lower', cmap=cm, vmin=lim[0], vmax=lim[1])
        cb = plt.colorbar(im, ax=ax[k], fraction=0.046, pad=0.03); cb.set_label(u)
        ax[k].set_title(f'Fig. 1E-{"abcd"[k]} – {n}'); ax[k].axis('off')
    fig.suptitle('Cartes ajustées (dti, non-linéaire Rician, options par défaut) – coupe axiale z = 24, sans overlay du masque', fontsize=9)
    plt.tight_layout(); plt.savefig(O + 'fig1E_maps.png', dpi=140); plt.close()
    fig, ax = plt.subplots(3, 3, figsize=(11, 11))
    for k in range(9):
        a = crop(np.where(m, F['D'][..., k], np.nan), m); lim = np.nanpercentile(np.abs(a), 99)
        im = ax.flat[k].imshow(a.T, origin='lower', cmap='RdBu_r', vmin=-lim, vmax=lim)
        cb = plt.colorbar(im, ax=ax.flat[k], fraction=0.046, pad=0.03); cb.set_label('µm²/ms')
        ax.flat[k].set_title(f'Fig. 1E-{"efghijklm"[k]} – D, volume {k + 1}'); ax.flat[k].axis('off')
    fig.suptitle('Les 9 volumes de la carte D (même coupe ; échelle bleu-blanc-rouge symétrique autour de 0 pour chaque volume)', fontsize=9)
    plt.tight_layout(); plt.savefig(O + 'fig1E_D.png', dpi=130); plt.close()
    rows = []
    for key, (i, j, z) in VOX.items():
        rows.append([key, F['FA'][i - 1, j - 1], *F['L'][i - 1, j - 1], *F['D'][i - 1, j - 1]])
    np.save(W + 'tab1F.npy', np.array(rows, dtype=object))
    for r in rows: print(r[0], ' '.join(f'{v:.3f}' for v in r[1:]))

def fig2A():
    names = ['Sphère (1, 1, 1)', 'Cigare (1.7, 0.3, 0.3)', 'Crêpe (1.5, 1.5, 0.3)']
    true_fa = [0, 0.799, 0.560]
    for k in range(1, 4):
        fig, ax = plt.subplots(figsize=(7.5, 4.3))
        L = sims[f'sh{k}_L']
        qplot(ax, sims[f'sh{k}_S'], sims[f'sh{k}_D'], 1.0,
              f'Fig. 2A-{"abc"[k-1]} – {names[k-1]}, SNR = 50 : FA vraie = {true_fa[k-1]:.3f}, FA ajustée = {L[3]:.3f}\n'
              f'L ajustées = {L[0]:.2f}, {L[1]:.2f}, {L[2]:.2f} µm²/ms   (o = signal simulé sans bruit, x = données bruitées, trait = fit)',
              true_S=sims[f'sh{k}_Sm'])
        plt.tight_layout(); plt.savefig(O + f'fig2A_{k}.png', dpi=140); plt.close()

def fig2B():
    for k, key in enumerate(VOX, 1):
        snr = int(sims[f'e{k}_snr'])
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
        L = sims[f'vox{k}_L']
        qplot(ax[0], sims[f'vox{k}_S'], sims[f'vox{k}_D'], sims[f'vox{k}_S0'], f'Données réelles – {NAMES[key]}\nFA = {L[3]:.2f} ; L = {L[0]:.2f}, {L[1]:.2f}, {L[2]:.2f}')
        L50 = sims[f'b{k}_50_L']
        qplot(ax[1], sims[f'b{k}_50_S'], sims[f'b{k}_50_D'], 1.0, f'Simulation, SNR = 50\nFA ajustée = {L50[3]:.2f} ; L = {L50[0]:.2f}, {L50[1]:.2f}, {L50[2]:.2f}', true_S=sims[f'b{k}_Sm'])
        Le = sims[f'e{k}_L']
        qplot(ax[2], sims[f'e{k}_S'], sims[f'e{k}_D'], 1.0, f'Simulation, SNR = {snr} (retenu)\nFA ajustée = {Le[3]:.2f} ; L = {Le[0]:.2f}, {Le[1]:.2f}, {Le[2]:.2f}', true_S=sims[f'b{k}_Sm'])
        fig.suptitle(f'Fig. 2B-{"abc"[k-1]} – {NAMES[key]} [{", ".join(map(str, VOX[key]))}] : simulation avec les L1, L2, L3 ajustés (1F) à côté du fit réel. '
                     f'Signal réel en unités arbitraires, simulé normalisé (S0 = 1).', fontsize=9)
        plt.tight_layout(); plt.savefig(O + f'fig2B_{key}.png', dpi=130); plt.close()

if __name__ == '__main__':
    which = sys.argv[1:] or ['1B', '1C', '1D', '2A', '2B']
    for w in which: globals()['fig' + w]()
