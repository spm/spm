function res = frequencies(this, varargin)
% Method for getting/setting frequencies of TF data
% FORMAT res = frequencies(this, varargin)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: frequencies.m 1270 2008-03-28 14:35:16Z stefan $

if isempty(varargin)
    if strcmp(transformtype(this), 'TF')
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
    sD = struct(this);
    sD.transform.ID = 'TF';
    sD.transform.frequencies = f;
    
    res = meeg(sD);
    
end
