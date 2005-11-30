function D = spm_eeg_inv_mesh_ui(S)

%=======================================================================
% Meshing user-interface routine
% commands the spatial normalization (if required) and the computation of
% the proper size individual size
%
% FORMAT D = spm_eeg_inv_mesh_ui(S)
% Input:
% S		    - input data struct (optional)
% Output:
% D			- same data struct including the meshing files and variables
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_mesh_ui.m 342 2005-11-30 18:27:52Z jeremie $

spm_defaults

if nargin == 0
    D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
elseif nargin == 1
    D = S;
else
	error(sprintf('Trouble reading the data file\n'));    
end

switch D.inv{end}.method
    
    case 'Imaging'
        val = length(D.inv);
        if val > 1
            D = spm_eeg_inv_copyfields(D,[1 0 0 0]);
        end
        
        Mflag = spm_input('Spatial Normalization?','+1','Yes|No',[1 2]);

        if Mflag == 1       % sMRI spatial normalization into MNI T1 template
            D = spm_eeg_inv_spatnorm(D);     
        end

        if Mflag == 2       % compute a template-based individual cortex mesh
            D = spm_eeg_inv_meshing(D);
            spm_eeg_inv_checkmeshes(D);
        end
        
    case 'ECD'
        val = length(D.inv);
        if val > 1
            D = spm_eeg_inv_copyfields(D,[1 0 0 0]);
        end
        [Pbin,Pinv_sn,kk,model] = spm_eeg_inv_model('Init',1);
        [pth,nam,ext,num] = spm_fileparts(model.fname);
        D.inv{end}.model = [nam,ext];
        [pth,namb,ext,num] = spm_fileparts(model.flags.fl_bin.Pimg);
        D.inv{end}.mesh.sMRI = [namb,ext];
        [pth,nam,ext,num] = spm_fileparts(Pbin(1,:));
        D.inv{end}.mesh.msk_iskusll = [nam,ext];
        [pth,nam,ext,num] = spm_fileparts(Pbin(2,:));
        D.inv{end}.mesh.msk_iscalp = [nam,ext];
        [pth,nam,ext,num] = spm_fileparts(Pinv_sn(1,:));
        D.inv{end}.mesh.invdef = [nam,ext];
        f_scVert = fullfile(pth,[namb,'_scVert','.mat']);
        scVert = model.head(end).XYZmm';
        save(f_scVert,'scVert')
        D.inv{end}.datareg.scalpvert = [namb,'_scVert','.mat'];
        
    otherwise
        error(sprintf('Trouble finding the method\n'));
end

save(D.fname,'D');