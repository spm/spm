function [U] = spm_get_ons(SPM,s)
% Returns input [designed effects] structures
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
%   U(i).P      - parameter struct.
%
%       U(i).P(p).name - parameter name
%       U(i).P(p).P    - parameter vector
%       U(i).P(p).h    - order of polynomial expansion
%       U(i).P(p).i    - sub-indices of u pertaining to P
%__________________________________________________________________________
%
% SLICE TIMING
%
% With longs TRs you may want to shift the regressors so that they are
% aligned to a particular slice. This is effected by resetting the
% values of defaults.stats.fmri.t and defaults.stats.fmri.t0 in
% spm_defaults. 
% defaults.stats.fmri.t is the number of time-bins per scan used when
% building regressors. Onsets are defined in temporal units of scans
% starting at 0.
% defaults.stats.fmri.t0 is the first time-bin at which the regressors are
% resampled to coincide with data acquisition. If defaults.stats.fmri.t0
% is set to 1 then the regressors will be appropriate for the first slice.
% If you want to temporally realign the regressors so that they match
% responses in the middle slice then make defaults.stats.fmri.t0 =
% defaults.stats.fmri.t/2 (assuming there is a negligible gap between
% volume acquisitions).
% Default values are defaults.stats.fmri.t=16 and defaults.stats.fmri.t0=1.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_get_ons.m 4083 2010-10-08 10:31:55Z guillaume $


% time units
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
        %------------------------------------------------------------------
        TR = T*dt;

    case 'secs'
        %------------------------------------------------------------------
        TR = 1;

    otherwise
        %------------------------------------------------------------------
        error('Unknown unit "%s".',UNITS);
end

% get inputs and names
%==========================================================================
try
    U   = SPM.Sess(s).U;
    v   = length(U);
catch

    %-prompt string
    %----------------------------------------------------------------------
    str = sprintf('Session %d: trial specification in %s',s,UNITS);
    spm_input(str,1,'d')

    U   = {};
    v   = spm_input('number of conditions/trials',2,'w1');
end

% get trials
%--------------------------------------------------------------------------
for i = 1:v

    % get names
    %----------------------------------------------------------------------
    try
        Uname     = U(i).name(1);
    catch
        str       = sprintf('name for condition/trial %d ?',i);
        Uname     = {spm_input(str,3,'s',sprintf('trial %d',i))};
        U(i).name = Uname;
    end

    % get main [trial] effects
    %======================================================================

    % onsets
    %----------------------------------------------------------------------
    try
        ons = U(i).ons;
        ons = ons(:);
    catch
        ons = [];
    end
    if isempty(ons)
        str      = ['vector of onsets - ' Uname{1}];
        ons      = spm_input(str,4,'r',' ',[Inf 1]);
        U(i).ons = ons(:);
    end

    % durations
    %----------------------------------------------------------------------
    try
        dur = U(i).dur;
        dur = dur(:);
    catch
        dur = [];
    end
    if isempty(dur)
        str = 'duration[s] (events = 0)';
        while 1
            dur = spm_input(str,5,'r',' ',[Inf 1]);
            if length(dur) == 1
                dur    = dur*ones(size(ons));
            end
            if length(dur) == length(ons), break, end
            str = sprintf('enter a scalar or [%d] vector',...
                length(ons));
        end
        U(i).dur = dur;
    end

    % peri-stimulus times {seconds}
    %----------------------------------------------------------------------
    pst   = [0:(k-1)]*T*dt - min(ons)*TR;
    for j = 1:length(ons)
        w      = [0:(k-1)]*T*dt - ons(j)*TR;
        v      = find(w >= 0);
        pst(v) = w(v);
    end


    % add parameters x trial interactions
    %======================================================================

    % get parameter stucture xP
    %----------------------------------------------------------------------
    try
        xP    = U(i).P;
        Pname = xP(1).name;

        switch Pname

            case 'none'
                %----------------------------------------------------------
                xP.name  = 'none';
                xP.h     = 0;

        end

    catch

        Pname = {'none','time','other'};
        Pname = spm_input('parametric modulation',6,'b',Pname);

        switch Pname

            case 'none'
                %----------------------------------------------------------
                xP(1).name  = 'none';
                xP(1).h     = 0;

            case 'time'   % units: minutes
                %----------------------------------------------------------
                xP(1).name  = 'time';
                xP(1).P     = ons*TR/60;
                xP(1).h     = spm_input('polynomial order',8,'n1',1);

            case 'other'
                %----------------------------------------------------------
                str   = ['# parameters (' Uname{1} ')'];
                for q = 1:spm_input(str,7,'n1',1);

                    % get names and parametric variates
                    %------------------------------------------------------
                    str   = sprintf('parameter %d name',q);
                    Pname = spm_input(str,7,'s');
                    P     = spm_input(Pname,7,'r',[],[length(ons),1]);

                    % order of polynomial expansion h
                    %------------------------------------------------------
                    h     = spm_input('polynomial order',8,'n1',1);

                    % sub-indices and inputs
                    %------------------------------------------------------
                    xP(q).name = Pname;
                    xP(q).P    = P(:);
                    xP(q).h    = h;

                end
        end % switch

    end % try

    % interaction with causes (u) - 1st = main effects
    %----------------------------------------------------------------------
    u     = ons.^0;
    for q = 1:length(xP)
        xP(q).i = [1, ([1:xP(q).h] + size(u,2))];
        for   j = 1:xP(q).h
            u   = [u xP(q).P.^j];
            str = sprintf('%sx%s^%d',Uname{1},xP(q).name,j);
            Uname{end + 1} = str;
        end
    end

    % orthogonalize inputs
    %----------------------------------------------------------------------
    u      = spm_orth(u);

    % and scale so sum(u*dt) = number of events, if event-related
    %----------------------------------------------------------------------
    if ~any(dur)
        u  = u/dt;
    end

    % create stimulus functions (32 bin offset)
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

    % place in ouputs structure
    %----------------------------------------------------------------------
    U(i).name = Uname;      % - input names
    U(i).dt   = dt;         % - time bin {seconds}
    U(i).u    = sf;         % - stimulus function matrix
    U(i).pst  = pst;        % - pst (seconds)
    U(i).P    = xP;         % - parameter struct

end % (v)
