%% GBM8360E Lab 3 - capture 1A (v4)
% Ouvre qMRLab 3.0.0 avec le modele dti et les donnees demo, puis enregistre
% la capture de la fenetre dans Downloads\fig1A_qMRLab.png
base = fullfile(getenv('USERPROFILE'),'Downloads');
q = fullfile(base,'qMRLab');
d = fullfile(base,'dti_demo');
if ~exist(fullfile(d,'DiffusionData.nii.gz'),'file'), unzip(fullfile(base,'dti.zip'), d); end
if ~exist(fullfile(d,'DiffusionData.nii'),'file'), gunzip(fullfile(d,'*.nii.gz'), d); end   % load_nii ne lit plus les .gz sous R2026
addpath(genpath(q));
% fermer les fenetres qMRLab precedentes et vider leur cache
delete(findall(groot,'Type','figure'));
s = getappdata(0); for f = fieldnames(s)', rmappdata(0, f{1}); end
% qMRLab 3.0.0 ne cree le navigateur de fichiers que pour le modele par defaut :
% on met dti comme modele par defaut (equivalent du bouton "Set Default")
Method = 'dti'; save(fullfile(q,'src','Common','Parameters','DefaultMethod.mat'),'Method');
data.DiffusionData = fullfile(d,'DiffusionData.nii');   % le GUI attend des chemins de fichiers
data.Mask          = fullfile(d,'Mask.nii');
app = qmrlab.gui.MainApp(dti, data);
drawnow; pause(10); drawnow;
exportapp(app.qMRILab, fullfile(base,'fig1A_qMRLab.png'));
disp('Capture enregistree : Downloads\fig1A_qMRLab.png');
