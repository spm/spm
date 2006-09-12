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
sdip = spm_eeg_inv_ecd_fitDip_ui(D);
D.inv{D.val}.inverse.sdip   = sdip;

% 2. if multiple seeds were used, summarise results by grouping solutions
resdip = spm_eeg_inv_ecd_sDipRes(sdip);
D.inv{D.val}.inverse.resdip = resdip;

% 3. display the results on the MRI
spm_eeg_inv_ecd_DrawDip('Init',resdip,D.inv{D.val}.mesh.sMRI)

save(D.fname,'D');
