F.mmc = '/dcm_fnirs/mmclab';
F.isomesh = '/dcm_fnirs/iso2mesh'; 
F.mesh = '/dcm_fnirs/MMC_Collins_Atlas_Mesh_Version_2L.mat'; 
F.atlas = '/dcm_fnirs/colin27_t1_tal_lin.nii'; 
F.sdpos = '/dcm_fnirs/pos_nirs.mat'; 

[G] = estimate_greens_mmclab(F); 
