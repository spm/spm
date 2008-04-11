function this = rmfield(this, stuff)
% Method for removing an object field
% FORMAT this = rmfield(this, stuff)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: rmfield.m 1373 2008-04-11 14:24:03Z spm $

eval(['obj = rmfield(obj.' stuff ');'])
