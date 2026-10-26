%% GBM8360E - Lab 3 - Exercices 1 et 2 (qMRLab 3.0.0, modele dti)
% Script d'appoint : il fait les calculs (masque, fits, valeurs aux voxels,
% simulations) et produit des figures de secours. Les captures d'ecran
% demandees dans l'enonce se prennent dans le GUI de qMRLab.
%
% A executer section par section (Ctrl+Entree) dans MATLAB R2026a/b.
% Toolboxes : Image Processing, Optimization, Statistics and ML.
%
% Eliott Dionnet - octobre 2026

%% 0. Reglages (A ADAPTER)
qMRLabDir = fullfile(getenv('HOME'),'qMRLab');      % dossier ou qMRLab est clone/dezippe
demoDir   = fullfile(pwd,'dti_demo','dti_data');     % dossier cree par "Download" dans le GUI
outDir    = fullfile(pwd,'lab3_out'); if ~exist(outDir,'dir'), mkdir(outDir); end

% Voxels choisis a la partie 1C [x y z] (coordonnees affichees dans le GUI) :
vox.genu       = [0 0 0];   % <-- A REMPLIR
vox.cortexGM   = [0 0 0];   % <-- A REMPLIR
vox.ventricule = [0 0 0];   % <-- A REMPLIR
sliceZ = vox.genu(3);       % coupe axiale utilisee en 1C/1D/1E/1F

addpath(genpath(qMRLabDir));

%% 1. Chargement des donnees demo + protocole
Model = dti;
% Si tu as sauvegarde le modele depuis le GUI (bouton Save), charge-le pour
% avoir exactement le meme protocole/options :
% Model = Model.loadObj(fullfile(demoDir,'dti.qmrlab.mat'));
data.DiffusionData = double(load_nii_data(fullfile(demoDir,'DiffusionData.nii.gz')));
nVol = size(data.DiffusionData,4);
assert(size(Model.Prot.DiffusionData.Mat,1)==nVol, ...
    'Le protocole (%d lignes) ne correspond pas aux donnees (%d volumes) : charger le protocole du GUI.', ...
    size(Model.Prot.DiffusionData.Mat,1), nVol);
Prot  = Model.Prot.DiffusionData.Mat;
Gnorm = Prot(:,4);
idx.b0    = find(Gnorm==0);
idx.b1000 = find(abs(Gnorm-0.0566)<1e-3);
idx.b1990 = find(abs(Gnorm-0.0800)<1e-3);
fprintf('Volumes : %d b=0, %d b=1000, %d b=1990 (sur %d)\n', numel(idx.b0), numel(idx.b1000), numel(idx.b1990), nVol);
fprintf('Taille des donnees : %s\n', mat2str(size(data.DiffusionData)));

%% 1B. Vues axiale / coronale / sagittale a b = 0, 1000, 1990 (figure de secours)
vols = [idx.b0(1) idx.b1000(1) idx.b1990(1)]; labs = {'b = 0','b = 1000','b = 1990'};
sz = size(data.DiffusionData); c = round(sz(1:3)/2);
figure('Name','1B','Color','w','Position',[50 50 1100 900]);
for k = 1:3
    V = data.DiffusionData(:,:,:,vols(k));
    lim = prctile(V(V>0),[1 99.5]);               % contraste ajuste par volume
    views = {rot90(V(:,:,c(3))), rot90(squeeze(V(:,c(2),:))), rot90(squeeze(V(c(1),:,:)))};
    names = {'axiale','coronale','sagittale'};
    for p = 1:3
        subplot(3,3,(k-1)*3+p); imagesc(views{p}, lim); axis image off; colormap gray;
        title(sprintf('Fig. 1B - %s, %s (vol %d)', names{p}, labs{k}, vols(k)));
    end
end
exportgraphics(gcf, fullfile(outDir,'fig1B_raw.png'), 'Resolution',150);

%% 1C. Fit voxel par voxel (memes resultats que "Fit voxel" du GUI)
names = fieldnames(vox);
for n = 1:numel(names)
    v = vox.(names{n});
    datavox = extractvoxel(data, v, {'DiffusionData'});
    FR = Model.fit(datavox);
    figure('Name',['1C ' names{n}],'Color','w'); Model.plotModel(FR, datavox);
    title(sprintf('Fig. 1C - %s [%d %d %d] : FA=%.2f, L=%.2f/%.2f/%.2f', names{n}, v, FR.FA, FR.L1, FR.L2, FR.L3));
    exportgraphics(gcf, fullfile(outDir,['fig1C_' names{n} '.png']), 'Resolution',150);
    fprintf('%-11s [%3d %3d %3d]  FA=%.3f  L1=%.3f L2=%.3f L3=%.3f um2/ms\n', names{n}, v, FR.FA, FR.L1, FR.L2, FR.L3);
end

%% 1D. Masque du cerveau sur UNE coupe axiale
% Methode : moyenne des b=0 -> seuil d'Otsu -> remplissage des trous ->
% plus grande composante connexe -> petite fermeture morphologique.
b0  = mean(data.DiffusionData(:,:,sliceZ,idx.b0), 4);
b0n = mat2gray(b0);
bw  = imbinarize(b0n, graythresh(b0n));
bw  = imfill(bw, 'holes');
bw  = bwareafilt(bw, 1);
bw  = imclose(bw, strel('disk',3));
bw  = imfill(bw, 'holes');
Mask = zeros(sz(1:3));  Mask(:,:,sliceZ) = bw;
save_nii_v2(Mask, fullfile(outDir,'Mask_slice.nii.gz'), fullfile(demoDir,'DiffusionData.nii.gz'), 2);
fprintf('Masque : %d voxels sur %d dans la coupe %d\n', nnz(bw), numel(bw), sliceZ);
figure('Name','1D','Color','w'); imagesc(rot90(b0)); axis image off; colormap gray; hold on;
visboundaries(rot90(bw),'Color','r'); title(sprintf('Fig. 1D - masque (contour rouge) sur b=0, coupe %d', sliceZ));
exportgraphics(gcf, fullfile(outDir,'fig1D_mask.png'), 'Resolution',150);
% -> Dans le GUI : charger lab3_out/Mask_slice.nii.gz dans le champ "Mask".

%% 1E. Fit de la coupe masquee + cartes
data.Mask = Mask;
FitResults = FitData(data, Model, 0);
save(fullfile(outDir,'FitResults_slice.mat'), '-struct', 'FitResults');
maps = {'FA',[0 1]; 'L1',[0 3]; 'L2',[0 3]; 'L3',[0 3]};
figure('Name','1E maps','Color','w','Position',[50 50 1200 330]);
for k = 1:4
    subplot(1,4,k); imagesc(rot90(FitResults.(maps{k,1})(:,:,sliceZ)), maps{k,2}); axis image off;
    colormap(gca, parula); cb = colorbar; if k>1, cb.Label.String = '\mum^2/ms'; end
    title(['Fig. 1E - ' maps{k,1}]);
end
exportgraphics(gcf, fullfile(outDir,'fig1E_maps.png'), 'Resolution',150);
figure('Name','1E D','Color','w','Position',[50 50 1100 900]);
Dlab = {'D(1)','D(2)','D(3)','D(4)','D(5)','D(6)','D(7)','D(8)','D(9)'};
for k = 1:9
    subplot(3,3,k); Dk = FitResults.D(:,:,sliceZ,k);
    m = max(abs(Dk(:))); imagesc(rot90(Dk), [-m m]); axis image off;
    colormap(gca, bwr_map()); colorbar; title(['Fig. 1E - ' Dlab{k}]);
end
exportgraphics(gcf, fullfile(outDir,'fig1E_D.png'), 'Resolution',150);

%% 1F. Valeurs ajustees aux 3 voxels
for n = 1:numel(names)
    v = vox.(names{n});
    Dv = squeeze(FitResults.D(v(1),v(2),v(3),:))';
    fprintf('\n%s [%d %d %d]\n  FA=%.3f  L1=%.3f  L2=%.3f  L3=%.3f\n  D = %s\n', names{n}, v, ...
        FitResults.FA(v(1),v(2),v(3)), FitResults.L1(v(1),v(2),v(3)), ...
        FitResults.L2(v(1),v(2),v(3)), FitResults.L3(v(1),v(2),v(3)), mat2str(Dv,3));
end

%% 2A. Simulations : sphere, cigare, crepe a SNR = 50
SimModel = dti;              % options et protocole par defaut
Opt.SNR  = 50;
shapes = {'Sphere',[1 1 1]; 'Cigare',[1.7 .3 .3]; 'Crepe',[1.5 1.5 .3]};
for s = 1:3
    L = shapes{s,2}; x = struct('L1',L(1),'L2',L(2),'L3',L(3));
    Lm = mean(L); FAtrue = sqrt(3/2)*sqrt(sum((L-Lm).^2))/sqrt(sum(L.^2));
    figure('Name',['2A ' shapes{s,1}],'Color','w');
    FR = SimModel.Sim_Single_Voxel_Curve(x, Opt);
    title(sprintf('Fig. 2A - %s : L = %s, FA vraie = %.3f, FA ajustee = %.3f', shapes{s,1}, mat2str(L), FAtrue, FR.FA));
    exportgraphics(gcf, fullfile(outDir,['fig2A_' shapes{s,1} '.png']), 'Resolution',150);
    fprintf('%-7s FA vraie = %.3f   FA ajustee = %.3f\n', shapes{s,1}, FAtrue, FR.FA);
end

%% 2B. Reproduire les voxels reels : simulation avec les L du 1F, plusieurs SNR
% Le protocole de simulation doit etre celui des donnees reelles (Model).
SNRs = [50 30 20 15 10];       % ajuster jusqu'a ce que la dispersion ressemble aux donnees
for n = 1:numel(names)
    v = vox.(names{n});
    x = struct('L1',FitResults.L1(v(1),v(2),v(3)),'L2',FitResults.L2(v(1),v(2),v(3)),'L3',FitResults.L3(v(1),v(2),v(3)));
    datavox = extractvoxel(data, v, {'DiffusionData'});
    FRreal = Model.fit(datavox);
    figure('Name',['2B ' names{n}],'Color','w','Position',[50 50 1500 400]);
    subplot(1,numel(SNRs)+1,1); Model.plotModel(FRreal, datavox); title(sprintf('Reel - %s', names{n}));
    for k = 1:numel(SNRs)
        subplot(1,numel(SNRs)+1,k+1); Opt.SNR = SNRs(k);
        Model.Sim_Single_Voxel_Curve(x, Opt); legend off; title(sprintf('Simulation SNR = %d', SNRs(k)));
    end
    exportgraphics(gcf, fullfile(outDir,['fig2B_' names{n} '.png']), 'Resolution',120);
end
% NB : la simulation normalise S0 a 1 alors que les donnees reelles sont en
% unites arbitraires : comparer la forme et la dispersion, pas l'echelle.

%% --- utilitaire : colormap divergente bleu-blanc-rouge
function cm = bwr_map(n)
if nargin<1, n = 256; end
h = n/2; t = linspace(0,1,h)';
cm = [[t t ones(h,1)]; [ones(h,1) flipud(t) flipud(t)]];
end
