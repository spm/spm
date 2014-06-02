function con = spm_dcm_contrasts(DCM,D)
% Make contrast vector for a DCM
% FORMAT con = spm_dcm_contrasts(DCM,D)
%
% DCM    - DCM structure or its filename
% D      - 'A','B' or 'C' i.e. connectivity matrix of interest
%
% con    - column vector specifying contrast weights
%__________________________________________________________________________
% Copyright (C) 2002-2014 Wellcome Trust Centre for Neuroimaging
 
% Will Penny
% $Id: spm_dcm_contrasts.m 6031 2014-06-02 12:49:52Z guillaume $
 

%-Load DCM if necessary
%--------------------------------------------------------------------------
if ~isstruct(DCM)
    load(DCM);
end
 
%-Ask user for contrast values
%--------------------------------------------------------------------------
con_struct  = spm_dcm_connectivity_ui(DCM,D,'Enter contrast for ');

%-Vectorize
%--------------------------------------------------------------------------
Ep          = DCM.Ep;      % MAP estimates
if ~isempty(con_struct)
    con     = spm_unvec(spm_vec(Ep)*0,Ep);
    con.(D) = con_struct.(D);
    con     = spm_vec(con);
end
