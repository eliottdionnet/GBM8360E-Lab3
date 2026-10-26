warning('off','all'); addpath(genpath('qMRLab')); addpath('lib');
M = dti; dd = double(load_nii_data('dti_demo/DiffusionData.nii.gz'));
V = [35 63 24; 28 78 24; 34 54 24];  % genou, cortex, ventricule
out = struct();
for k=1:3
  v=V(k,:); s.DiffusionData = squeeze(dd(v(1),v(2),v(3),:)); r = dti_fit_oct(M,s,true);
  sig = computesigmanoise(M.Prot.DiffusionData.Mat, s.DiffusionData);
  out.(sprintf('vox%d_S',k)) = s.DiffusionData; out.(sprintf('vox%d_D',k)) = r.D;
  out.(sprintf('vox%d_L',k)) = [r.L1 r.L2 r.L3 r.FA]; out.(sprintf('vox%d_sigma',k)) = sig; out.(sprintf('vox%d_S0',k)) = r.S0;
  printf('vox %d [%d %d %d] FA %.3f L %.3f %.3f %.3f  S0 %.0f sigma %.1f SNR %.1f\n', k, v, r.FA, r.L1, r.L2, r.L3, r.S0, sig, r.S0/sig);
end
% 2A
randn('seed',42); rand('seed',42);
shapes = [1 1 1; 1.7 .3 .3; 1.5 1.5 .3];
for k=1:3
  x.L1=shapes(k,1); x.L2=shapes(k,2); x.L3=shapes(k,3); Sm = M.equation(x); sg = max(Sm)/50;
  d.DiffusionData = ricernd(Sm, sg); r = dti_fit_oct(M,d,true);
  out.(sprintf('sh%d_Sm',k))=Sm; out.(sprintf('sh%d_S',k))=d.DiffusionData; out.(sprintf('sh%d_D',k))=r.D; out.(sprintf('sh%d_L',k))=[r.L1 r.L2 r.L3 r.FA];
  printf('shape %d FA fit %.3f L %.3f %.3f %.3f\n', k, r.FA, r.L1, r.L2, r.L3);
end
% 2B
SNRs = [50 30 20 15 10];
for k=1:3
  Lv = out.(sprintf('vox%d_L',k)); x.L1=Lv(1); x.L2=Lv(2); x.L3=Lv(3); Sm = M.equation(x);
  for n=1:numel(SNRs)
    randn('seed',100+10*k+n); rand('seed',100+10*k+n);
    d.DiffusionData = ricernd(Sm, max(Sm)/SNRs(n)); r = dti_fit_oct(M,d,true);
    out.(sprintf('b%d_%d_S',k,SNRs(n)))=d.DiffusionData; out.(sprintf('b%d_%d_D',k,SNRs(n)))=r.D; out.(sprintf('b%d_%d_L',k,SNRs(n)))=[r.L1 r.L2 r.L3 r.FA];
    printf('2B vox %d SNR %d: fit FA %.3f L %.3f %.3f %.3f\n', k, SNRs(n), r.FA, r.L1, r.L2, r.L3);
  end
  out.(sprintf('b%d_Sm',k)) = Sm;
  snr = round(out.(sprintf('vox%d_S0',k))/out.(sprintf('vox%d_sigma',k)));
  randn('seed',200+k); rand('seed',200+k);
  d.DiffusionData = ricernd(Sm, max(Sm)/snr); r = dti_fit_oct(M,d,true);
  out.(sprintf('e%d_S',k))=d.DiffusionData; out.(sprintf('e%d_D',k))=r.D; out.(sprintf('e%d_L',k))=[r.L1 r.L2 r.L3 r.FA]; out.(sprintf('e%d_snr',k))=snr;
  printf('2B vox %d SNR estime %d: fit FA %.3f L %.3f %.3f %.3f\n', k, snr, r.FA, r.L1, r.L2, r.L3);
end
save('-v7','out/sims.mat','-struct','out');
