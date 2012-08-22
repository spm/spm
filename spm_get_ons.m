function [U] = spm_get_ons(SPM,s)
% Return input [designed effects] structure
% FORMAT [U] = spm_get_ons(SPM,s)
%
% SPM   - SPM structure (see spm_fMRI_design.m)
% s     - session number
%
% U     - (1 x n)   struct array of (n) trial-specific structures
%
%   U(i).name   - cell of names for each input or cause
%   U(i).u      - inputs or stimulus function matrix
%   U(i).dt     - time bin (seconds)
%   U(i).ons    - onsets    (in SPM.xBF.UNITS)
%   U(i).dur    - durations (in SPM.xBF.UNITS)
%   U(i).orth   - orthogonalise inputs?
%   U(i).P      - parameter structure
%
%       U(i).P(p).name - parameter name
%       U(i).P(p).P    - parameter vector
%       U(i).P(p).h    - order of polynomial expansion
%       U(i).P(p).i    - sub-indices of u pertaining to P
%__________________________________________________________________________
%
% Note on Slice Timing:
%
% With longs TRs you may want to shift the regressors so that they are
% aligned to a particular slice. This is controlled by two variables:
% fMRI_T is the number of time-bins per scan used when building regressors.
% Onsets are defined in temporal units of scans starting at 0.
% fMRI_T0 is the first time-bin at which the regressors are resampled to 
% coincide with data acquisition. If fMRI_T0 is set to 1 then the 
% regressors will be appropriate for the first slice.
% If you want to temporally realign the regressors so that they match
% responses in the middle slice then make fMRI_T0 = fMRI_T/2 (assuming
% there is a negligible gap between volume acquisitions).
% Default values are defined in spm_defaults.m
%__________________________________________________________________________
% Copyright (C) 1999-2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_get_ons.m 4855 2012-08-22 15:50:22Z guillaume $


%-Time units
%--------------------------------------------------------------------------
k     = SPM.nscan(s);
T     = SPM.xBF.T;
dt    = SPM.xBF.dt;
try
    UNITS = SPM.xBF.UNITS;
catch
    UNITS = 'scans';
end
switch UNITS
    case 'scans'
        TR = T*dt;
    case 'secs'
        TR = 1;
    otherwise
        error('Unknown unit "%s".',UNITS);
end

%-Get inputs and names
%==========================================================================
U   = SPM.Sess(s).U;

%-Get trials
%--------------------------------------------------------------------------
for i = 1:numel(U)

    %-Get main [trial] effects
    %======================================================================

    %-Names
    Uname = U(i).name(1);
    
    %-Onsets
    %----------------------------------------------------------------------
    ons   = U(i).ons(:);

    %-Durations
    %----------------------------------------------------------------------
    dur   = U(i).dur(:);

    %-Peri-stimulus times {seconds}
    %----------------------------------------------------------------------
    pst   = [0:(k-1)]*T*dt - min(ons)*TR;
    for j = 1:length(ons)
        w      = [0:(k-1)]*T*dt - ons(j)*TR;
        v      = find(w >= 0);
        pst(v) = w(v);
    end

    %-Add parameters x trial interactions
    %======================================================================

    %-Get parameter structure xP
    %----------------------------------------------------------------------
    xP    = U(i).P;
    if strcmpi(xP(1).name,'none'), xP.h = 0; end

    %-Interaction with causes (u) - 1st = main effects
    %----------------------------------------------------------------------
    u     = ons.^0;
    for q = 1:length(xP)
        xP(q).i = [1, ([1:xP(q).h] + size(u,2))];
        for   j = 1:xP(q).h
            u   = [u xP(q).P.^j];
            Uname{end + 1} = sprintf('%sx%s^%d',Uname{1},xP(q).name,j);
        end
    end

    %-Orthogonalise inputs
    %----------------------------------------------------------------------
    if ~isfield(U(i),'orth') || U(i).orth
        u     = spm_orth(u);
    end

    %-And scale so sum(u*dt) = number of events, if event-related
    %----------------------------------------------------------------------
    if ~any(dur)
        u     = u/dt;
    end

    %-Create stimulus functions (32 bin offset)
    %======================================================================
    ton       = round(ons*TR/dt) + 33;               % onsets
    tof       = round(dur*TR/dt) + ton + 1;          % offset
    sf        = sparse((k*T + 128),size(u,2));
    ton       = max(ton,1);
    tof       = max(tof,1);
    for j = 1:length(ton)
        if size(sf,1) > ton(j)
            sf(ton(j),:) = sf(ton(j),:) + u(j,:);
        end
        if size(sf,1) > tof(j)
            sf(tof(j),:) = sf(tof(j),:) - u(j,:);
        end
    end
    sf        = cumsum(sf);                         % integrate
    sf        = sf(1:(k*T + 32),:);                 % stimulus

    %-Place in ouputs structure U
    %----------------------------------------------------------------------
    U(i).name = Uname;      % - input names
    U(i).dt   = dt;         % - time bin {seconds}
    U(i).u    = sf;         % - stimulus function matrix
    U(i).pst  = pst;        % - pst (seconds)
    U(i).P    = xP;         % - parameter struct

end
