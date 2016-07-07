function res = sensors(this, type, newsens)
% Sets and gets sensor fields for EEG and MEG
% returns empty matrix if no sensors are defined.
% FORMAT res = sensors(this, type, newsens)
%   type - 'EEG' or 'MEG'
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: sensors.m 6829 2016-07-07 10:16:46Z vladimir $

if nargin<2
    error('Sensor type (EEG or MEG) must be specified');
end

switch lower(type)
    case 'eeg'
        if nargin < 3
            if isfield(this.sensors, 'eeg')
                res = this.sensors.eeg;
            else
                res = [];
            end
        else
            this.sensors(1).eeg = newsens;
            res = check(this);
        end
    case 'meg'
        if nargin < 3
            if isfield(this.sensors, 'meg')
                res = this.sensors.meg;
            else
                res = [];
            end
        else
            this.sensors(1).meg = newsens;
            res = check(this);
        end
    case 'src'
        if nargin < 3
            if isfield(this.sensors, 'src')
                res = this.sensors.src;
            else
                res = [];
            end
        else
            this.sensors(1).src = newsens;
            res = check(this);
        end
    otherwise
        error('Unsupported sensor type');
end

if  nargin < 3 && ~isempty(res) && ismember(lower(type), {'eeg', 'meg'}) && this.montage.Mind > 0
    sens = res;
    montage = this.montage.M(this.montage.Mind);
    if  ~isempty(intersect(sens.label, montage.labelorg))
        sensmontage = montage;
        [sel1, sel2] = spm_match_str(sens.label, sensmontage.labelorg);
        sensmontage.labelorg = sensmontage.labelorg(sel2);
        sensmontage.tra = sensmontage.tra(:, sel2);
        selempty  = find(all(sensmontage.tra == 0, 2));
        sensmontage.tra(selempty, :) = [];
        sensmontage.labelnew(selempty) = [];
        
        if isfield(sensmontage, 'chantypeorg')
            sensmontage.chantypeorg = sensmontage.chantypeorg(sel2);
        end
        if isfield(sensmontage, 'chanunitorg')
            sensmontage.chanunitorg = sensmontage.chanunitorg(sel2);
        end
        if isfield(sensmontage, 'chantypenew')
            sensmontage.chantypenew(selempty) = [];
        end
        if isfield(sensmontage, 'chanunitnew')
            sensmontage.chanunitnew(selempty) = [];
        end
        
        chanunitorig = sens.chanunit(sel1);
        chantypeorig = sens.chantype(sel1);
        labelorg     = sens.label;
        
        keepunused   = 'no'; % not sure if this is good for all cases
        
        sens = ft_apply_montage(sens, sensmontage, 'keepunused', keepunused);
        
        if strcmpi(type, 'MEG')
            if isfield(sens, 'balance') && ~isequal(sens.balance.current, 'none')
                balance = ft_apply_montage(getfield(sens.balance, sens.balance.current), sensmontage, 'keepunused', keepunused);
            else
                balance = sensmontage;
            end
            
            sens.balance.custom = balance;
            sens.balance.current = 'custom';
        end
        
        
        % If all the original channels contributing to a new channel have
        % the same units, transfer them to the new channel. This might be
        % wrong if the montage itself changes the units by scaling the data.
        if isfield(sens, 'chanunit')
            chanunit = sens.chanunit;
        else
            chanunit = repmat({'unknown'}, numel(sens.label), 1);
        end
        
        if isfield(sens, 'chantype')
            chantype = sens.chantype;
        else
            chantype = repmat({'unknown'}, numel(sens.label), 1);
        end
        
        for j = 1:numel(sens.label)
            k = strmatch(sens.label{j}, sensmontage.labelnew, 'exact');
            if ~isempty(k)
                if isequal(chanunit{j}, 'unknown')
                    unit = unique(chanunitorig(~~abs(sensmontage.tra(k, :))));
                    if numel(unit)==1
                        chanunit(j) = unit;
                    elseif strcmpi(type, 'MEG')
                        chanunit{j} = 'T';
                    else
                        chanunit{j} = 'V';
                    end
                end
                
                if isequal(chantype{j}, 'unknown')
                    ctype = unique(chantypeorig(~~abs(sensmontage.tra(k, :))));
                    if numel(ctype)==1
                        chantype(j) = ctype;
                    elseif strcmpi(type, 'MEG')
                        chantype{j} = 'megmag';
                    else
                        chantype{j} = 'eeg';
                    end
                end
            else %channel was not in the montage, but just copied
                k = strmatch(sens.label{j}, labelorg, 'exact');
                if isequal(chanunit{j}, 'unknown')
                    chanunit(j) = chanunitorig(k);
                end
                if isequal(chantype{j}, 'unknown')
                    chantype(j) = chantypeorig(k);
                end
            end
        end
        
        sens.chanunit = chanunit;
        sens.chantype = chantype;
    end
    
    if strcmpi(type, 'MEG')
        res = ft_datatype_sens(sens, 'amplitude', 'T', 'distance', 'mm');
    else
        res = ft_datatype_sens(sens, 'amplitude', 'V', 'distance', 'mm');
    end
end
