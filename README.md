# GBM8360E – Lab 3 (A2026) – Eliott Dionnet

Code et données du laboratoire 3 (qMRLab 3.0.0, modèle `dti`), Polytechnique Montréal.

## Contenu

| Dossier / fichier | Rôle |
|---|---|
| `code/lab3_ex1_ex2.m` | Exercices 1 et 2 sous MATLAB : masque, fits voxel et coupe, valeurs aux voxels, simulations |
| `code/lab3_capture_1A.m` | Ouvre le GUI qMRLab avec le modèle dti et les données démo (capture 1A) |
| `code/ex1_ex2_octave/` | Version utilisée pour les résultats du rapport (code qMRLab exécuté sous Octave) + figures Python ; voir son README |
| `code/lab3_ex3.m` | Exercice 3 : reproduction du problème, diagnostic, réparation (2 méthodes) |
| `code/repair_ex3.py` | Exercice 3 : réparation en Python (nibabel) |
| `code/figures_ex3.py`, `code/fit_maps_octave.m`, `code/sim_water_octave.m` | Figures de l'exercice 3 |
| `data/Mask_slice.nii.gz` | Masque d'une coupe (exercice 1D, coupe z = 24) |
| `data/collab_dwi_fixed.nii.gz` | Exercice 3 : données du collaborateur avec les volumes remis dans l'ordre (méthode 1) |
| `data/collab_dti_reversedprot.qmrlab.mat`, `data/collab_dwi_reversedprot.bval/.bvec` | Exercice 3 : protocole remis dans l'ordre des données (méthode 2) |

## Exercice 3 en une phrase

Les 19 volumes de `collab_dwi.nii.gz` étaient stockés dans l'ordre inverse du protocole (le vrai b = 0 était le dernier volume) ;
inverser l'axe des volumes (ou l'ordre du protocole) donne L ≈ 2 µm²/ms et FA ≈ 0.06 dans l'eau.

## Dépendances

qMRLab 3.0.0 (https://github.com/qMRLab/qMRLab), MATLAB R2026b (Image Processing, Optimization, Statistics) ou Octave 8.4 ;
Python 3 avec numpy, scipy, nibabel, matplotlib, scikit-image.

Les données démo (`dti.zip`) se téléchargent depuis qMRLab (bouton « Download example ») ; les données de l'exercice 3 depuis
https://github.com/GBM8360/Lab3_data/releases/tag/1.0.0.

Réalisé avec l'aide de Claude (Anthropic).
