function [res] = montage(this,action,varargin)
% Method for specifying an online montage, or setting one to use
% FORMAT
%   res = montage(this, 'add', montage)
%           Adding a montage to the meeg object, see format here under
%
%   res = montage(this, 'action', idx)
%           Setting, checking, getting or removing a montage in the object,
%           depending on the action string and index idx of montage.
% Actions:
% - add        -> adding a montage to the object
% - switch     -> switch between montage, 0 being no applied montage 
%                 (switch to 0 by default if no index passed)
% - remove     -> removing montage, one at a time or any list.
% - getnumber  -> returning the number of montage(s) available
% - getindex   -> return current montage index
% - getname    -> returning a list of montage name (by default the current
%                 one if no list is passed)
% - getmontage -> returning the current or any other montage structure, 
%                 depending on list provided (current one by default if 
%                 no list passed).
% _______________________________________________________________________
% Copyright (C) 2011-2017 Wellcome Trust Centre for Neuroimaging

% Remy Lehembre & Christophe Phillips
% Cyclotron Research Centre, University of Liege, Belgium
% $Id: montage.m 7111 2017-06-16 09:01:09Z guillaume $

% Montage definition in the object structure by simply adding a 'montage'
% field in the object structure:
% D.montage.M(1:k)
%          .Mind   - 0       => no montage applied
%                  \ 1...N   => use any of the N specified montage
% where M is a structure with
%  - name      [optional]
%  - tra       M*N montage matrix
%  - labelnew  M*1 cell array of labels
%  - labelorg  N*1 cell array of labels
%  - channels  [optional] same format as the main 'channels' field in the
%                           meeg object
%      * label
%      * bad
%      * type
%      * x_plot2D
%      * y_plot2D
%      * units
%
% Note: If no channels information is passed, then we'll try to guess what
% to put in. It's easy for simple channel selection and re-referencing but
% not possible for more general cases.

switch lower(action)
    
    case 'add'
        % adding a montage to the object
        % structure is passed, add montage
        mont = varargin{1};
        %check that all info are consistent
        if size(mont.tra,1)~=length(mont.labelnew) || ...
                size(mont.tra,2)~=length(mont.labelorg)
            error('Montage Matrix inconsistent with original or new number of electrodes.')
        end
        if size(mont.labelorg,1)~=size(mont.tra,2)
            mont.labelorg = mont.labelorg';
        end
        if size(mont.labelnew,1)~=size(mont.tra,1)
            mont.labelnew = mont.labelnew';
        end
        
        % Check if there are already some montages, if not initialize
        % => set "ind" as the i^th montage
        this = struct(this);
        if ~isfield(this,'montage')
            this.montage = [];
            this.montage.M = [];
            ind = 1;
        else
            if isfield(this.montage,'M')
                ind = length(this.montage.M)+1;
            else
                this.montage = [];
                this.montage.M = [];
                ind = 1;
            end
        end
        
        % write montage info to fields M and Mind of montage
        if isfield(mont,'name') && ~isempty(mont.name)
            this.montage.M(ind).name = mont.name;
        else
            this.montage.M(ind).name = ['montage #',num2str(ind)];
        end
        this.montage.M(ind).tra = mont.tra;
        this.montage.M(ind).labelnew = mont.labelnew;
        this.montage.M(ind).labelorg = mont.labelorg;
        this.montage.Mind = ind;
        
        % fill channel information
        if isfield(mont,'channels')
            % use channel information provided
            % NO check performed here !!!
            this.montage.M(ind).channels = mont.channels;
        else
            % try to derive it from object channel info
            disp('No new channels information : setting channels info automatically.')
            this = set_channels(this,mont);
        end
        res = meeg(this);
        
    case 'switch'
        % switching between different montages
        if nargin==2
            idx = 0; % by default, no montage applied.
        else
            idx = varargin{1};
        end
        
        if ischar(idx) && ~isempty(this.montage.M)
            idx = strmatch(idx, {this.montage.M.name});
            if isempty(idx), idx = -1; end                
        end
        
        if idx>numel(this.montage.M) || idx<0
            error('Specified montage index is erroneous.')
        else
            this.montage.Mind = idx;
        end
        res = meeg(this);
        
    case {'remove', 'clear'} 
        if strcmpi(action, 'clear')
            idx = 1:numel(this.montage.M);
        else
            % removing one or more montages
            if nargin==2
                idx = this.montage.Mind; % by default, removing current montage.
            else
                idx = varargin{1};
            end
        end
        
        if ~isempty(idx)
            if any(idx>numel(this.montage.M) | idx<0)
                error('Specified montage index is erroneous.')
            else
                this.montage.M(idx) = [];
                if any(idx==this.montage.Mind)
                    % removing current montage -> no montage applied
                    this.montage.Mind = 0;
                elseif any(idx < this.montage.Mind)
                    % removing another montage, keep current (adjusted) index
                    this.montage.Mind = this.montage.Mind - ...
                        sum(idx < this.montage.Mind) ;
                end
            end
        end
        res = meeg(this);
        
    case 'getnumber'
        % getting the total numbr of available montages
        res = numel(this.montage.M);
        
    case 'getindex'
        % getting the index of the current montage selected
        res = this.montage.Mind;
        
    case 'getname'
        % getting the name of current or any other montage(s)
        if nargin==2
            idx = this.montage.Mind; % by default, current montage.
        else
            idx = varargin{1};
        end
        
        if ~isnumeric(idx)
            idx = 1:numel(this.montage.M);
        end
        
        if any(idx>numel(this.montage.M)) || any(idx<0)
            error('Specified montage index is erroneous.')
        elseif isequal(idx, 0) || isempty(this.montage.M)
            res = 'none';
        else
            res = {this.montage.M(idx).name};
            if numel(res) == 1
                res = char(res);
            end
        end
        
    case 'getmontage'
        % getting the current montage structure or any other one
        if nargin==2
            idx = this.montage.Mind; % by default, current montage.
        else
            idx = varargin{1};
        end
        if idx>numel(this.montage.M) || idx<0
            error('Specified montage index is erroneous.')
        elseif idx==0
            res = [];
        else
            res = this.montage.M(idx);
        end
    otherwise
        % nothing planned for this case...
        error('Wrong use of the ''montage'' method.')
end


%==========================================================================
function S = set_channels(S,mont)

% provided montage does not necessarily apply to all channels
% this adjusts it to include also the unused channels
mont_orig = [];
mont_orig.labelorg = {S.channels(:).label};
mont_orig.labelnew = mont_orig.labelorg;
mont_orig.tra = eye(numel(mont_orig.labelnew));
mont         = ft_apply_montage(mont_orig, mont);

% set channels "default" value, try to guess values from main channels
% definition in the object.

% Use new channel labels and set bad=0
idx = S.montage.Mind;
for ii=1:length(mont.labelnew)
    S.montage.M(idx).channels(ii).label = mont.labelnew{ii};
    S.montage.M(idx).channels(ii).bad = 0;
end

% Set new electrodes as bad if they include a bad channel
res = [S.channels.bad];
res = find(res);

newbads = find(any(mont.tra(:, res),2));
for ii=1:length(newbads)
    S.montage.M(idx).channels(newbads(ii)).bad = 1;
end

% set channel info: type, units
l_EEG = [];
for ii=1:length(S.montage.M(idx).channels)
    l_chan_org = find(mont.tra(ii,:));
    % 'type'
    type_ii = unique({S.channels(l_chan_org).type});
    if numel(type_ii)==1
        S.montage.M(idx).channels(ii).type = type_ii{1};
    else
        % mixing different types of channels
        S.montage.M(idx).channels(ii).type = 'Other';
    end
    l_EEG = [l_EEG ii]; % list EEG channels
    
    % 'units'
    units_ii = unique({S.channels(l_chan_org).units});
    if numel(units_ii)==1
        S.montage.M(idx).channels(ii).units = units_ii{1};
    else
        % mixing different units of channels
        S.montage.M(idx).channels(ii).units = 'unknown';
    end
end

% Deal with "new" channel positions:
% - when there is one channel with weight >0, assume it's a simple channel
%   selection or re-referencing
%   => use the info from channel with weight>0
% - otherwise reset channel position to []

for ii=1:length(S.montage.M(idx).channels)
    tra_ii = mont.tra(ii,:);
    
    if sum(tra_ii>0)==1
        % 1 channel extracted or re-ref -> keep coord of (+) channel
        S.montage.M(idx).channels(ii).X_plot2D = ...
            S.channels(tra_ii>0).X_plot2D;
        S.montage.M(idx).channels(ii).Y_plot2D = ...
            S.channels(tra_ii>0).Y_plot2D;
    else
        % more than 1 channel with >0 weight, don't know what to do...
        S.montage.M(idx).channels(ii).X_plot2D = [];
        S.montage.M(idx).channels(ii).Y_plot2D = [];
    end
end    
