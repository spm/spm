function res = frequencies(this, varargin)
% Method for getting/setting frequencies of TF data
% FORMAT res = frequencies(this, varargin)
% _________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: frequencies.m 3196 2009-06-11 12:54:47Z vladimir $

if isempty(varargin)
    if strncmpi(transformtype(this), 'TF',2)
        res = getset(this, 'transform', 'frequencies');
    else
        res = [];
    end
else
    f = varargin{1};
        
    if any(f) <= 0 || any(~isnumeric(f))
        error('Frequencies must be positive numbers'); res = []; return
    end

    % can't use getset because both information must be set at the same
    % time
    this.transform.ID = 'TF'; % could it be TFphase?
    this.transform.frequencies = f;
    
    res = this;
    
end
