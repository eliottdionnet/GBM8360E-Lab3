warning('off','all'); addpath(genpath('qMRLab')); 
D='lab3-for_students/';
M = dti; M = M.loadObj([D 'collab_dti.qmrlab.mat']);
rand('seed',1); randn('seed',1);
x.L1=2; x.L2=2; x.L3=2; Sm = M.equation(x); sig = max(Sm)/50; dat.DiffusionData = ricernd(Sm,sig); r = M.fit(dat);
Ssim = dat.DiffusionData; Dsim = r.D; Lsim=[r.L1 r.L2 r.L3 r.FA];
save('-v7','octave_out/sim_water.mat','Sm','Ssim','Dsim','Lsim');
disp(Lsim)
