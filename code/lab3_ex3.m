%% GBM8360E - Lab 3 - Exercice 3 : fantome du collaborateur
% Reproduit le probleme, teste l'hypothese (volumes stockes dans l'ordre
% inverse du protocole) et repare le jeu de donnees de deux facons.
% Ecrit pour MATLAB R2026 + qMRLab 3.0.0 (calculs verifies sous Octave 8.4).

%% 0. Reglages
qMRLabDir = fullfile(getenv('HOME'),'qMRLab');
D   = fullfile(pwd,'lab3-for_students');          % dossier dezippe
out = fullfile(pwd,'lab3_ex3_out'); if ~exist(out,'dir'), mkdir(out); end
addpath(genpath(qMRLabDir));

Model = dti;
Model = Model.loadObj(fullfile(D,'collab_dti.qmrlab.mat'));   % protocole + options du collaborateur
data.DiffusionData = double(load_nii_data(fullfile(D,'collab_dwi.nii.gz')));
data.Mask          = double(load_nii_data(fullfile(D,'collab_mask.nii.gz')));
Prot = Model.Prot.DiffusionData.Mat;

% Voxels de controle [x y z] (coordonnees MATLAB, 1-based)
vox = struct('eau',[12 20 1], 'tube_x',[16 32 1], 'tube_y',[32 48 1], 'croisement',[32 32 1]);

%% A. Reproduire le probleme
Fit0 = FitData(data, Model, 0);
afficher(Fit0, vox, 'Donnees recues');

%% B. Indices
% B1 : intensite moyenne de chaque volume dans le fantome vs etiquette du protocole
m = data.Mask(:,:,1) > 0;
moy = squeeze(sum(sum(data.DiffusionData(:,:,1,:).*m,1),2)) / nnz(m);
bval = round(Prot(:,4).^2 / 0.0566^2 * 100) * 10;           % b approx. a partir de Gnorm
disp(table((1:size(Prot,1))', bval, Prot(:,1:3), moy, 'VariableNames',{'vol','b_etiquete','dir','signal_moyen'}));
% -> le volume etiquete b=0 n'est pas le plus brillant ; le dernier (etiquete b=1500) l'est.

% B2 : simulation d'un voxel d'eau attendu (L1=L2=L3=2) avec le protocole du collaborateur
Opt.SNR = 50; x = struct('L1',2,'L2',2,'L3',2);
figure('Color','w'); Model.Sim_Single_Voxel_Curve(x, Opt); title('Simulation : eau, D = 2 \mum^2/ms, SNR 50');
figure('Color','w'); dv = extractvoxel(data, vox.eau, {'DiffusionData'}); Model.plotModel(Model.fit(dv), dv);
title('Voxel d''eau reel (donnees recues)');

%% C/D. Reparation, methode 1 : remettre les volumes dans l'ordre du protocole
fixed = data; fixed.DiffusionData = data.DiffusionData(:,:,:,end:-1:1);
save_nii_v2(fixed.DiffusionData, fullfile(out,'collab_dwi_fixed.nii.gz'), fullfile(D,'collab_dwi.nii.gz'), 64);
Fit1 = FitData(fixed, Model, 0);
afficher(Fit1, vox, 'Methode 1 : volumes inverses');

%% D. Reparation, methode 2 : garder les donnees, inverser l'ordre des lignes du protocole
Model2 = Model; Model2.Prot.DiffusionData.Mat = Prot(end:-1:1,:); Model2 = Model2.UpdateFields;
Model2.saveObj(fullfile(out,'collab_dti_reversedprot.qmrlab.mat'));
Fit2 = FitData(data, Model2, 0);
afficher(Fit2, vox, 'Methode 2 : protocole inverse');
fprintf('Ecart max methode 1 vs 2 (FA) : %g\n', max(abs(Fit1.FA(:)-Fit2.FA(:))));

%% E. Cartes avant / apres
figure('Color','w','Position',[50 50 1300 600]);
cartes = {'FA','L1','L2','L3'};
for k = 1:4
    subplot(2,4,k);   imagesc(rot90(Fit0.(cartes{k})), [-3 3]); axis image off; colorbar; colormap(gca, bwr());
    title([cartes{k} ' - recu']);
    subplot(2,4,k+4); imagesc(rot90(Fit1.(cartes{k})), [0 2.5]*(k>1) + [0 1]*(k==1)); axis image off; colorbar; colormap(gca, parula);
    title([cartes{k} ' - repare']);
end

%% ---------------------------------------------------------------------
function afficher(F, vox, titre)
fprintf('\n== %s ==\n', titre);
n = fieldnames(vox);
for k = 1:numel(n)
    v = vox.(n{k});
    fprintf('%-11s [%2d %2d]  L = %6.2f %6.2f %6.2f um2/ms   FA = %.2f\n', n{k}, v(1), v(2), ...
        F.L1(v(1),v(2),v(3)), F.L2(v(1),v(2),v(3)), F.L3(v(1),v(2),v(3)), F.FA(v(1),v(2),v(3)));
end
end

function cm = bwr(n)
if nargin<1, n = 256; end
h = n/2; t = linspace(0,1,h)';
cm = [[t t ones(h,1)]; [ones(h,1) flipud(t) flipud(t)]];
end
