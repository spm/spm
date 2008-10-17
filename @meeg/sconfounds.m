function res = sconfounds(this, newsconfounds)
% Method for getting/setting spatial confounds
% FORMAT res = sconfounds(this, newsconfounds)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: sconfounds.m 2355 2008-10-17 13:38:02Z vladimir $

if nargin == 2
    [sel1, sel2] = spm_match_str(chanlabels(this), newsconfounds.label);

    if length(sel1)~=numel(meegchannels(this))
        error('The spatial confounds do not match the MEEG channels.');
    end
    
    if any(newsconfounds.bad(sel2)) && ~all(badchannels(this, sel1(find(newsconfounds.bad(sel2)))))
        warning('Setting additional channels to bad to match spatial confounds.');
        this = badchannels(this, find(newsconfounds.bad(sel2)), 1);
    end

    newsconfounds.label = newsconfounds.label(sel2);
    newsconfounds.coeff = newsconfounds.coeff(sel2, :);
    newsconfounds.bad = newsconfounds.bad(sel2);

    this.other.sconfounds = newsconfounds;
    
    res = this;
else
    chanind = setdiff(meegchannels(this), badchannels(this));
    if ~isfield(this, 'sconfounds')
        res = zeros(length(chanind), 1);
        return;
    end    
    
    [sel1, sel2] = spm_match_str(chanlabels(this, chanind), this.other.sconfounds.label);
    
    if any(this.other.sconfounds.bad(sel2))
        error(['Channels ' sprintf('%s ', this.other.sconfounds.label{sel2}) ' should be set to bad.']);
    end
    
    res = this.other.sconfounds.coeff(sel2, :);
end
