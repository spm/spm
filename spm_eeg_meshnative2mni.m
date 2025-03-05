function nativemesh=spm_eeg_meshmni2native(mnimesh,mesh)
%function nativemesh=spm_eeg_meshmni2native(mnimesh,mesh)
%% Uses mesh ( spm mesh structure D.inv{:}.mesh ) to compute transform to 
%% express mnimesh in native space
%% replicates code segment from headmodel section of SPM code


defs.comp{1}.inv.comp{1}.def = {mesh.def};
defs.comp{1}.inv.space = {mesh.sMRI};
defs.out{1}.surf.surface = {mnimesh};
defs.out{1}.surf.savedir.savesrc = 1;
out = spm_deformations(defs);
nativemesh=gifti(out.surf{1});
