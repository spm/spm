function D = spm_eeg_inv_ecd_ui(S)
% Fits dipole(s) onto a bit of EEG data.
% To do all this, I use bits of a routine previously written for the DipFit
% toolbox. I try to render things a bit more coherent...
%
% FORMAT D = spm_eeg_inv_ecd_ui(S)
% Input:
% S		    - input data struct (optional)
% Output:
% D			- same data struct including the details of the model
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips
% $Id$

% 1. call dipfit gui to get all the parameters and fit the dipoles
% 2. if multiple seeds were used, summarise the results by grouping the
% solutions
% 3. display the results on the MRI

spm_defaults

if nargin == 0
    D = spm_eeg_ldata;
elseif nargin == 1
    D = S;
else
	error(sprintf('Trouble reading the data file\n'));    
end

% 1. call dipfit gui to get all the parameters and fit the dipoles
[sdip,fit_opt,Psave] = spm_eeg_inv_ecd_fitDip_ui(D);
if ~isfield(D.inv{end}.inverse,'sdip') || isempty(D.inv{end}.inverse.sdip)
    D.inv{end}.inverse.sdip{1} = Psave;
else
    D.inv{end}.inverse.sdip{end+1} = Psave;
end

% 2. if multiple seeds were used, summarise the results by grouping the
% solutions
[resdip,Pres] = spm_eeg_inv_ecd_sDipRes(sdip);
if ~isfield(D.inv{end}.inverse,'resdip') || isempty(D.inv{end}.inverse.resdip)
    D.inv{end}.inverse.resdip{1} = Pres;
else
    D.inv{end}.inverse.resdip{end+1} = Pres;
end

% 3. display the results on the MRI
spm_eeg_inv_ecd_DrawDip('Init',resdip,D.inv{end}.mesh.sMRI)

save(D.fname,'D');
