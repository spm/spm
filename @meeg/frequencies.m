function res = frequencies(this, ind, f)
% Method for getting/setting frequencies of TF data
% FORMAT res = frequencies(this, ind, values)
% _________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: frequencies.m 5079 2012-11-25 18:38:18Z vladimir $

if nargin >1        
    if ~isnumeric(ind)
        ind = 1:nfrequencies(this);
    end
end

if nargin < 3
    if strncmpi(transformtype(this), 'TF',2)
        res = this.transform.frequencies;
    else
        res = [];
        return
    end
    if exist('ind', 'var') == 1
        res = res(ind);
    end    
else       
     if ~strncmpi(transformtype(this), 'TF',2)
         error('Frequencies can only be assigned to a TF dataset');
     end
    
    if any(f) <= 0 || any(~isnumeric(f))
        error('Frequencies must be positive numbers');
    end
        
    if length(ind)~=length(f) || max(ind)>size(this, 2)
          error('Wrong frequency axis or indices'); 
    end

    if length(ind) == size(this.data, 2)
        this.transform.frequencies = f;
    else
        this.transform.frequencies(ind) = f;
    end
    
    res = this;    
end
