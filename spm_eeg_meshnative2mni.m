function mnimesh=spm_eeg_meshnative2mni(nativemesh,mesh)
%%function mnimesh=spm_eeg_meshnative2mni(nativemesh,mesh)
%% Uses mesh ( spm mesh structure D.inv{:}.mesh ) to compute transform to 
%% express
%% nativemesh(gifti filename) in native MRI space
%% as
%% mnimesh as gitfi structure in mni space
%% replicates code segment from headmodel section of SPM code


defs.comp{1}.inv.comp{1}.def = {mesh.def};
defs.comp{1}.inv.space = {mesh.sMRI};
defs.out{1}.surf.surface = {nativemesh};
defs.out{1}.surf.savedir.savesrc = 1;
out = spm_deformations(defs);
mnimesh=gifti(out.surf{1});
