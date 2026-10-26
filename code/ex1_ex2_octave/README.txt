Exercices 1 et 2 - calculs faits avec le code qMRLab 3.0.0 sous Octave 8.4
(dossiers attendus : qMRLab/ (clone v3.0.0), dti_demo/ (zip demo dezippe), out/, figs/)

1. prot.txt : protocole par defaut du modele dti (save -ascii de dti().Prot.DiffusionData.Mat)
2. python mask.py           -> out/Mask_slice.nii.gz (coupe z = 24)
3. octave ex1fit.m          -> out/fit_slice.mat (fit non-lineaire Rician de la coupe masquee)
4. octave ex12sim.m         -> out/sims.mat (fits des 3 voxels, simulations 2A et 2B)
5. python figs12.py         -> figures 1B, 1C, 1D, 1E, 2A, 2B

lib/dti_fit_oct.m : copie de dti.fit ou l'etape non-lineaire tourne aussi sous Octave
lib/pdf.m         : pdf('rician',...) (absent d'Octave sans le package statistics)
Version MATLAB equivalente : ../lab3_ex1_ex2.m
