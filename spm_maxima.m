function spm_maxima(SPM,VOL,hReg)
% List of local maxima for a single suprathreshold cluster
% FORMAT spm_maxima(SPM,VOL,hReg)
%
% SPM    - structure containing SPM, distribution & filtering details
%        - See spm_list('ListCluster',... for list of required fields.
%
% VOL    - structure containing details of volume analysed
%        - See spm_list('ListCluster',... for list of required fields.
%
% hReg   - handle of MIP XYZ registry object (see spm_XYZreg for details)
%_______________________________________________________________________
%
% spm_maxima functionality has been absorbed into spm_list.
%
% See spm_list('ListCluster',SPM,VOL)
%_______________________________________________________________________
% %W% Karl Friston, Andrew Holmes %E%


%-Print warning of obsolescence
%-----------------------------------------------------------------------
warning('spm_maxima is grandfathered, use spm_list(''ListCluster'',... instead')


%-Pass arguments on to spm_list
%-----------------------------------------------------------------------
spm_list('ListCluster',SPM,VOL);
