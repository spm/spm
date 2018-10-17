function labels = spm_create_labels(S)
% Create n  numbered labels using a base string as template
% FORMAT labels = spm_create_labels(S)
%   S           - input structure
%   Fields of S:
%   S.base      - Template string     - Default: 'T'
%   S.n         - number of labels    - Default:  1 
%
% Output:
%  labels       - cell array of labels
%
% Example:
%       S = [];
%       S.base = 'TRIG';
%       S.n = 100;
%       labels = spm_create_labels(S);
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_create_labels.m 7452 2018-10-17 15:29:54Z tim $

%-Set default values
%--------------------------------------------------------------------------
if ~nargin, S = struct(); end
if ~isfield(S, 'base'), S.base = 'T'; end
if ~isfield(S, 'n'),    S.n    = 1;   end

%-Create labels
%--------------------------------------------------------------------------
pad = numel(num2str(S.n));
fmt = [S.base,'%0',num2str(pad),'d'];
labels = arrayfun(@(x) sprintf(fmt,x), 1:S.n, 'UniformOutput',false);
