function str = spm_file(str,varargin)
% Character array (or cell array of strings) handling facility
% FORMAT str = spm_file(str,option)
% str        - character array, or cell array of strings
% option     - string of requested item - one among: {'path', 'cpath', 
%              'fpath', 'basename', 'ext', 'filename', 'number'}
%
% FORMAT str = spm_file(str,opt_key,opt_val,...)
% str        - character array, or cell array of strings
% opt_key    - string of targeted item - one among {'path', 'basename',
%              'ext', 'filename', 'number', 'prefix', 'suffix'}
% opt_val    - string of new value for feature
%__________________________________________________________________________
%
% Definition of items:
%
% <cpath> = <fpath>filesep<filename>
% <filename> = <basename>.<ext><number>
% <path> = empty or full path or relative path
%__________________________________________________________________________
%
% Examples:
%
% spm_file('C:\data\myimage.nii', 'prefix','rp_', 'ext','.txt')
% returns 'C:\data\rp_myimage.txt' on a Windows platform
%
% spm_file({'/home/karl/software/spm8/spm.m'},'path','/home/karl/spm12')
% returns {'/home/karl/spm12/spm.m'}
%
% spm_file('/home/karl/software/spm8/spm.m','filename')
% returns 'spm.m', and
% spm_file('/home/karl/software/spm8/spm.m','basename')
% returns 'spm'
%
% spm_file('SPM.mat','fpath')
% returns '/home/karl/data/stats' (i.e. pwd), while
% spm_file('SPM.mat','path')
% returns '', and
% spm_file('SPM.mat','cpath')
% returns '/home/karl/data/stats/SPM.mat'
%__________________________________________________________________________
%
% See also: spm_fileparts, spm_select, spm_file_ext, spm_existfile
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_file.m 4445 2011-08-26 17:53:00Z guillaume $


needchar = ischar(str);
options = varargin;

str = cellstr(str);

%-Get item
%==========================================================================
if numel(options) == 1
    for n=1:numel(str)
        [pth,nam,ext,num] = spm_fileparts(str{n});
        switch lower(options{1})
            case 'path'
                str{n} = pth;
            case 'basename'
                str{n} = nam;
            case 'ext'
                str{n} = ext(2:end);
            case 'filename'
                str{n} = [nam ext num];
            case 'number'
                str{n} = num;
            case 'cpath'
                str{n} = spm_select('CPath',str{n});
            case 'fpath'
                str{n} = spm_fileparts(spm_select('CPath',str{n}));
            otherwise
                error('Unknown option.');
        end
    end
    options = {};
end

%-Set item
%==========================================================================
while ~isempty(options)
    for n=1:numel(str)
        [pth,nam,ext,num] = spm_fileparts(str{n});
        switch lower(options{1})
            case 'path'
                pth = options{2};
            case 'basename'
                nam = options{2};
            case 'ext'
                ext = options{2};
                if ~isempty(ext) && ext(1) ~= '.'
                    ext = ['.' ext];
                end
                num = '';
            case 'filename'
                nam = options{2};
                ext = '';
            case 'number'
                num = options{2};
            case 'prefix'
                nam = [options{2} nam];
            case 'suffix'
                nam = [nam options{2}];
            otherwise
                warning('Unknown item ''%s'': ignored',lower(options{1}));
        end
        str{n} = fullfile(pth,[nam ext num]);
    end
    options([1 2]) = [];
end

if needchar
    str = char(str);
end
