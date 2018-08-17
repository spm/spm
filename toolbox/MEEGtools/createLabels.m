function labels = createLabels(S)
% Create n  numbered labels using a base string as template
% FORMAT labels = createLabels(S)
%   S           - input structure
%   Fields of S:
%   S.base      - Template string     - Default: 'T'
%   S.n         - number of labels    - Default:  1 
%
% Output:
%  labels      - MEEG object (also written to disk)
%
% Example:
%       S =[];
%       S.base = 'TRIG';
%       S.n =100;
%       labels = createLabels(S);
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id$


%-Set default values
%--------------------------------------------------------------------------
if ~isfield(S, 'base'),        S.base   = 'T'; end
if ~isfield(S, 'n'),           S.n   = 1; end

%- Create labels
%--------------------------------------------------------------------------

    pad = numel(num2str(S.n));
    cmd = [S.base,'%0',num2str(pad),'d'];
    labels = {};
    for i = 1:S.n
        labels{i} = sprintf(cmd,i);
    end
end