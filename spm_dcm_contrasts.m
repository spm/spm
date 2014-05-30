function [con] = spm_dcm_contrasts(DCM_filename,D)
% Make contrast vector for a DCM model
% FORMAT [con] = spm_dcm_contrasts(DCM_filename,D)
%
% DCM_filename  DCM file name
% D             'A','B' or 'C' i.e. connectivity matrix of interest
%
% con           Column vector specifying contrast weights
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Will Penny
% $Id: spm_dcm_contrasts.m 6026 2014-05-30 11:09:03Z peter $
 
% Load DCM if necessary
%--------------------------------------------------------------------------
try
    P   = DCM_filename;
    load(P);
catch
    DCM = DCM_filename;
end
 
% Ask user for contrast values
%--------------------------------------------------------------------------
con_struct = spm_dcm_connectivity_ui(DCM,D,'Enter contrast for ');

% Vectorize
%--------------------------------------------------------------------------
Ep      = DCM.Ep;           % MAP estimates
if ~isempty(con_struct)
    con     = spm_unvec(spm_vec(Ep)*0,Ep);
    con.(D) = con_struct.(D);
    con     = spm_vec(con);
end
