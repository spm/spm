function res = sconfounds(this, newsconfounds, append)
% Method for getting/setting spatial confounds
% FORMAT res = sconfounds(this, newsconfounds)
% _______________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: sconfounds.m 6437 2015-05-14 12:27:21Z vladimir $

if nargin >= 2
    meegind = indchantype(this, 'MEEG');
    
    [sel1, sel2] = match_str(chanlabels(this, meegind), newsconfounds.label);

    sel1 = meegind(sel1);        
    
    if length(sel1)<length(meegind)
        error('The spatial confounds do not match the MEEG channels.');
    end
    
    if any(newsconfounds.bad(sel2)) && ~all(badchannels(this, sel1(find(newsconfounds.bad(sel2)))))
        warning('Setting additional channels to bad to match spatial confounds.');
        this = badchannels(this, find(newsconfounds.bad(sel2)), 1);
    end

    newsconfounds.label = newsconfounds.label(sel2);
    newsconfounds.coeff = newsconfounds.coeff(sel2, :);
    newsconfounds.bad = newsconfounds.bad(sel2);
    
    if nargin == 3 && isfield(this.other, 'sconfounds')
        oldsconfounds = this.other.sconfounds;
        
        if ~isequal(newsconfounds.label, oldsconfounds.label)
            error('New confounds incompatible with old. Cannot append.');
        end
        
        newsconfounds.coeff = [oldsconfounds.coeff newsconfounds.coeff];
    end            

    this.other(1).sconfounds = newsconfounds;
    
    res = this;
else
    chanind = indchantype(this, 'MEEG', 'GOOD');
    if ~isfield(this, 'sconfounds')
        res = zeros(length(chanind), 1);
        return;
    end    
    
    [sel1, sel2] = match_str(chanlabels(this, chanind), this.other.sconfounds.label);
    
    if any(this.other.sconfounds.bad(sel2))
        error(['Channels ' sprintf('%s ', this.other.sconfounds.label{sel2}) ' should be set to bad.']);
    end
    
    res = this.other.sconfounds.coeff(sel2, :);
end
