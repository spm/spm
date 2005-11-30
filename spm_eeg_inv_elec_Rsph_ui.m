function D = spm_eeg_inv_elec_Rsph_ui(S)

% Deals with the solution of the forward problem:
% There are 2 seperate cases:
% 1. individual head model, based on the subject anatomy
% 2. standard head model, based on the template
% 
% For an individual head model :
% - project the sensor location, on the scalp surface (EEG), not sure what
% to do for MEG...
% - Estimate the best fitting spheres, and prepare the realistic sphere
% approach
% For a standard head model :
% - load the head model built on the template in MNI space (actually the single subject
% canonical image)
% -select your standard electrode setup
%
% To do all this, I use bits of a routine previously written for the DipFit
% toolbox. I try to render things a bit more coherent...
%
% FORMAT D = spm_eeg_inv_elec_Rsph_ui(S)
% Input:
% S		    - input data struct (optional)
% Output:
% D			- same data struct including the meshing files and variables
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips
% $Id$

spm_defaults

if nargin == 0
    D = spm_eeg_ldata;
elseif nargin == 1
    D = S;
else
	error(sprintf('Trouble reading the data file\n'));    
end

switch D.inv{end}.method
    case 'Imaging'
        % Not sure Jeremie will use this routine but...
        % then it would be easy to make all the rest available to 'imaging'
        % source reconstruction.
        
    case 'ECD'
        % Process data for ECD
        q_indiv = spm_input('Model type :','+1','individual|standard',[1,2],1)
        
        if q_indiv==1
            D = elec_Rsph(D);
        else
            D = std_model(D);
        end
    otherwise
        error(sprintf('Trouble finding the method\n'));
end



%_________________________________________________________________________
%_________________________________________________________________________
function D = elec_Rsph(D);
% Deal with the individual head model

% Start loading the various bits: model and electrodes location
load(D.inv{end}.model)
load(D.inv{end}.datareg.sens_coreg)
% Then put the electrodes in the model
flags_el = struct('q_RealLoc',1,'br_only',model.flags.fl_tess.br_only);
[electrodes,flags_el] = ...
    spm_eeg_inv_model('Elec2Scalp',model.head(end),sensreg',[],flags_el);
model.electrodes = electrodes ;
model.flags.fl_elec = flags_el ;
% How define the name, and so the order, of the channels ???
% Let's assume the order of the channels in the polhemus file is 
% the same as the order of the acquisition of the data
% => use the data structure to give channels names.
if length(D.channels.eeg)==model.electrodes.nr
    model.electrodes.names = D.channels.name(D.channels.eeg)
else
    spm('alert!',strvcat('Sorry different # of channels in data and sensors locations you provided',...
        'I''ll proceed but be aware there is a mismatch!'),'Channels problem');
end

% Set the sphere model
flags_Rs.br_only = model.flags.fl_tess.br_only;
[spheres,a,b,c,flags_Rs] = spm_eeg_inv_Rsph(model,[],flags_Rs);
model.spheres = spheres;
model.flags.fl_RealS = flags_Rs;
% Save things in the end
save(D.inv{end}.model,'model')

%_________________________________________________________________________
%_________________________________________________________________________
function D = std_model(D)
% Deal with the standard head model

% Put the model into the M/EEG data structure
load('standard_head_model.mat')
[pth,nam,ext] = fileparts(D.fname);
Pmod = fullfile(pth,[nam,'_std_model.mat']);
model.fname = Pmod
D.inv{end}.model = model.fname;
D.inv{end}.mesh.sMRI = model.flags.fl_bin.Pimg;
D.inv{end}.mesh.invdef = model.flags.fl_bin.Pinv_sn(1,:);
f_scVert = fullfile(pth,[nam,'_scVert','.mat']);
scVert = model.head(end).XYZmm';
save(f_scVert,'scVert')
D.inv{end}.datareg.scalpvert = f_scVert;

% Input standard electrode sets
[set_Nel,set_name] = spm_eeg_inv_electrset;
el_set = spm_input('Which set of electrodes', ...
					'+1','m',set_name);
[el_sphc,el_name] = spm_eeg_inv_electrset(el_set) ;
flags_el.q_RealLoc = 0;
flags_el.q_RealNI = 0;
[electrodes,flags_el] = spm_eeg_inv_model('Elec2Scalp',model.head(end), ...
                                    el_sphc,el_name,flags_el);
model.electrodes = electrodes;
model.flags.fl_elec = flags_el;

% Set the sphere model
flags_Rs.br_only = 0;
[spheres,a,b,c,flags_Rs] = spm_eeg_inv_Rsph(model,[],flags_Rs);
model.spheres = spheres;
model.flags.fl_RealS = flags_Rs;
% Save things in the end
save(D.inv{end}.model,'model')
        
        