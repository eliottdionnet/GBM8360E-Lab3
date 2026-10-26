warning('off','all'); addpath(genpath('qMRLab')); 
D='lab3-for_students/';
Model = dti; Model = Model.loadObj([D 'collab_dti.qmrlab.mat']);
dd = double(load_nii_data([D 'collab_dwi.nii.gz'])); mk=double(load_nii_data([D 'collab_mask.nii.gz']));
names={'orig','fixed'};
for h=1:2
  if h==2, dd = dd(:,:,:,end:-1:1); end
  L=zeros(64,64,3); FA=zeros(64,64); Dm=zeros(64,64,9);
  for i=1:64, for j=1:64
    if mk(i,j)
      vox.DiffusionData = squeeze(dd(i,j,1,:)); r=Model.fit(vox);
      L(i,j,:)=[r.L1 r.L2 r.L3]; FA(i,j)=r.FA; Dm(i,j,:)=r.D;
    end
  end, end
  save('-v7',['octave_out/map_' names{h} '.mat'],'L','FA','Dm');
end
