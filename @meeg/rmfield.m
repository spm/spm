function this = rmfield(this, stuff)
% Method for removing an object field
% FORMAT this = rmfield(this, stuff)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

eval(['obj = rmfield(obj.' stuff ');'])
