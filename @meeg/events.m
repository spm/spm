function res = events(this, varargin)
% Method for getting/setting events per trial
% FORMAT res = events(this, ind, event)
%   ind = indices of trials
% _______________________________________________________________________
% Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: events.m 5613 2013-08-15 11:56:07Z vladimir $

if nargin == 2
    res = getset(this, 'trials', 'events', varargin{:});
elseif nargin == 3 && ischar(varargin{2})
    ev     = getset(this, 'trials', 'events', varargin{1});
    onsets = trialonset(this, varargin{1});
    if ~iscell(ev)
        ev = {ev};
        strct = 1;
    else
        strct = 0;
    end
    
    for j = 1:numel(ev)
        event = ev{j};
        if ~isempty(event)
            for i = 1:numel(event)
                event(i).time = event(i).time - onsets(j);
                
                if isequal(varargin{2}, 'samples')
                    if onsets(j) == 0
                        event(i).sample = event(i).time*this.Fsample;
                    else
                        event(i).sample = event(i).time*this.Fsample+1;
                    end
                    
                    event(i).duration = ceil(event(i).duration*this.Fsample);
                    
                    event(i).sample   = max(round(event(i).sample), 1);
                end
            end
            if isequal(varargin{2}, 'samples')
                event = rmfield(event, 'time');
            end
        end
        ev{j} = event;
    end
    
    if strct
        res = ev{1};
    else
        res = ev;
    end
else
    res = getset(this, 'trials', 'events', varargin{:});
end