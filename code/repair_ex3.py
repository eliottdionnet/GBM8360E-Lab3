"""GBM8360E Lab 3 - Exercice 3 : reparation du jeu de donnees du collaborateur.

Les 19 volumes de collab_dwi.nii.gz sont stockes dans l'ordre INVERSE du
protocole (bval/bvec et collab_dti.qmrlab.mat). Deux reparations equivalentes :
  1) inverser l'axe des volumes du NIfTI (protocole inchange)      -> collab_dwi_fixed.nii.gz
  2) garder le NIfTI et inverser l'ordre des bval/bvec               -> *_reversedprot.bval/.bvec
     (pour qMRLab : inverser les lignes de Prot, cf. lab3_ex3.m)

Usage : python repair_ex3.py <dossier lab3-for_students> <dossier de sortie>
"""
import sys, os
import numpy as np
import nibabel as nib

src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)

img = nib.load(os.path.join(src, 'collab_dwi.nii.gz'))
data = np.asarray(img.dataobj)
assert data.shape[-1] == 19, data.shape

# Methode 1 : volumes remis dans l'ordre du protocole (meme en-tete, meme type)
fixed = nib.Nifti1Image(data[..., ::-1].copy(), img.affine, img.header)
nib.save(fixed, os.path.join(out, 'collab_dwi_fixed.nii.gz'))

# Methode 2 : protocole remis dans l'ordre des donnees
bval = np.loadtxt(os.path.join(src, 'collab_dwi.bval'))
bvec = np.loadtxt(os.path.join(src, 'collab_dwi.bvec'))
np.savetxt(os.path.join(out, 'collab_dwi_reversedprot.bval'), bval[None, ::-1], fmt='%d')
np.savetxt(os.path.join(out, 'collab_dwi_reversedprot.bvec'), bvec[:, ::-1], fmt='%.6f')

# Controle rapide : dans le fantome, le volume b=0 doit etre le plus brillant
mask = np.asarray(nib.load(os.path.join(src, 'collab_mask.nii.gz')).dataobj).reshape(data.shape[:2]) > 0
m_before = data[:, :, 0, :][mask].mean(0)
m_after = data[:, :, 0, ::-1][mask].mean(0)
print('volume le plus brillant avant :', m_before.argmax() + 1, '(etiquete b =', int(bval[m_before.argmax()]), ')')
print('volume le plus brillant apres :', m_after.argmax() + 1, '(etiquete b =', int(bval[m_after.argmax()]), ')')
