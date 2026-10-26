import nibabel as nib, numpy as np
from skimage.filters import threshold_otsu
from scipy import ndimage as ndi
D='dti_demo/'; O='out/'
Z=23  # 0-based  (MATLAB z = 24)
img=nib.load(D+'DiffusionData.nii.gz')
G=np.array([0.0]*0)
P=np.loadtxt('out/prot.txt')
b0idx=np.where(P[:,3]==0)[0]
b0=np.asarray(img.dataobj[:,:,Z,:])[...,b0idx].mean(-1).astype(float)
t=threshold_otsu(b0); bw=b0>t
bw=ndi.binary_fill_holes(bw)
lab,n=ndi.label(bw); sizes=ndi.sum(bw,lab,range(1,n+1)); bw=lab==(np.argmax(sizes)+1)
st=np.zeros((7,7),bool); yy,xx=np.mgrid[-3:4,-3:4]; st=(xx**2+yy**2)<=9
bw=ndi.binary_closing(bw,st); bw=ndi.binary_fill_holes(bw)
M=np.zeros(img.shape[:3],np.uint8); M[:,:,Z]=bw
h=img.header.copy(); h.set_data_dtype(np.uint8); h.set_data_shape(M.shape)
nib.save(nib.Nifti1Image(M,img.affine,h),O+'Mask_slice.nii.gz')
np.save(O+'b0_slice.npy',b0)
print('otsu',t,'mask voxels',bw.sum(),'of',bw.size)
