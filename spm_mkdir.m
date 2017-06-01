function sts = spm_mkdir(varargin)
% Make new directory trees
% FORMAT sts = spm_mkdir(dir,...)
% dir    - character array, or cell array of strings
%
% sts    - true if all directories were successfully created or already
%          existing, false otherwise.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mkdir.m 7088 2017-06-01 16:18:44Z guillaume $


sts = true;
if nargin > 0
    d1 = cellstr(varargin{1});
    for i=1:numel(d1)
        if ~exist(spm_select('cpath',d1{i}),'dir')
            status = mkdir(d1{i});
            sts = sts & status;
        end
        if nargin > 1
            d2 = cellstr(varargin{2});
            for j=1:numel(d2)
                status = spm_mkdir(fullfile(d1{i},d2{j}),varargin{3:end});
                sts = sts & status;
            end
        end
    end
else
    error('Not enough input arguments.');
end
