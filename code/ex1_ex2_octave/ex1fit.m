warning('off','all'); addpath(genpath('qMRLab')); addpath('lib');
M = dti; z = 24;
dd = double(load_nii_data('dti_demo/DiffusionData.nii.gz'));
mk = load_nii_data('out/Mask_slice.nii.gz');
[nx,ny,~,nv] = size(dd);
FA=zeros(nx,ny); L=zeros(nx,ny,3); D=zeros(nx,ny,9); S0=zeros(nx,ny); tic;
for i=1:nx, for j=1:ny
  if mk(i,j,z)
    s.DiffusionData = squeeze(dd(i,j,z,:)); r = dti_fit_oct(M,s,true);
    FA(i,j)=r.FA; L(i,j,:)=[r.L1 r.L2 r.L3]; D(i,j,:)=r.D; S0(i,j)=r.S0;
  end
end, end
toc
save('-v7','out/fit_slice.mat','FA','L','D','S0');
