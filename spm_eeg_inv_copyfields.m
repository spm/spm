function D = spm_eeg_inv_copyfields(S,Cflags,Vcheck)

%=======================================================================
% Copy needed fields from the first inverse analysis on the same data
%
% FORMAT D = spm_eeg_inv_copyfields(S,Cflags,Vcheck)
% Input:
% S		    - input data struct (optional)
% Cflags    - indicators of fields to duplicate
%             (default: minimum files for further analysis on the same
%             subject)
% Vcheck    - check the emptiness of the strucutre to be filled in
%             1: checking    (not empty -> no copy)
%             0: no checking (copy anyway)
%
% Output:
% D			- same data struct including the new files and parameters
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_copyfields.m 621 2006-09-12 17:22:42Z karl $

def_Cflags = [1 1 1 0];
def_Vcheck = 1;

if nargin == 0
    D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
    D = spm_eeg_ldata(D);
    Cflags = def_Cflags;
    Vcheck = def_Vcheck;
elseif nargin == 1
    if isstruct(S);
        D = S;
        Cflags = def_Cflags;
        Vcheck = def_Vcheck;
    else
        error(sprintf('Wrong input arguments\n'));
    end
    Cflags = def_Cflags;
elseif nargin == 2
    D = S;
    Vcheck = def_Vcheck;
elseif nargin == 3;
    D = S;
else
    error(sprintf('Error: wrong input arguments\n'));
end

if ~isfield(D,'inv')
    error(sprintf('Error: no inverse structure for these data\n'));
end

try
    val = D.val;
catch
    val = length(D.inv);
end

% Meshing
if Cflags(1) ~= 0
    names_m = [];
    names_m = fieldnames(D.inv{val}.mesh);
    count_m = 0;
    if Vcheck == 1
        % check for empty structure
        for i = 1:length(names_m)
            count_m = count_m + isempty(getfield(D.inv{val}.mesh,names_m{i}));
        end
        if count_m ~= length(names_m)
            disp('Mesh structure is not empty - no copy has been made for those fields');
        end
    end
    if Vcheck == 0 | count_m == length(names_m);
        % copy fields according to options
        if Cflags(1) == 2
            D.inv{val}.mesh = D.inv{val-1}.mesh;
        elseif Cflags(1) == 1
            D.inv{val}.mesh.sMRI        = D.inv{val-1}.mesh.sMRI;
            D.inv{val}.mesh.nobias      = D.inv{val-1}.mesh.nobias;
            D.inv{val}.mesh.def         = D.inv{val-1}.mesh.def;
            D.inv{val}.mesh.invdef      = D.inv{val-1}.mesh.invdef;
            D.inv{val}.mesh.msk_iskull  = D.inv{val-1}.mesh.msk_iskull;
            D.inv{val}.mesh.msk_scalp   = D.inv{val-1}.mesh.msk_scalp;
        end
        disp('Mesh fields copied');
    end
end

% Data registration
if Cflags(2) ~= 0
    names_d = [];
    names_d = fieldnames(D.inv{val}.datareg);
    count_d = 0;
    if Vcheck == 1
        % check for empty structure
        for i = 1:length(names_d)
            count_d = count_d + isempty(getfield(D.inv{val}.datareg,names_d{i}));
        end
        if count_d ~= length(names_d)
            disp('Datareg structure is not empty - no copy has been made for those fields');
        end
    end
    if Vcheck == 0 | count_d == length(names_d)
        % copy fields according to options
        if Cflags(2) == 2
            D.inv{val}.datareg = D.inv{val-1}.datareg;
        elseif Cflags(2) == 1
            D.inv{val}.datareg.sens     = D.inv{val-1}.datareg.sens;
        end
        disp('Datareg fields copied');
    end
end

% Forward parameters
if Cflags(3) ~= 0
    names_f = [];
    names_f = fieldnames(D.inv{val}.forward);
    count_f = 0;
    if Vcheck == 1
        % check for empty structure
        for i = 1:length(names_f)
            count_f = count_f + isempty(getfield(D.inv{val}.forward,names_f{i}));
        end
        if count_f ~= length(names_f)
            disp('Forward structure is not empty - no copy has been made for those fields');
        end
    end
    if Vcheck == 0 | count_f == length(names_f)
        % copy fields according to options
        if Cflags(3) == 2
            D.inv{val}.forward = D.inv{val-1}.forward;
        elseif Cflags(2) == 1
            D.inv{val}.forward.bst_channel     = D.inv{val-1}.datareg.bst_channel;
        end
        disp('Forward fields copied');
    end
end

% Inverse parameters
if Cflags(4) ~= 0
    names_i = [];
    names_i = fieldnames(D.inv{val}.inverse);
    count_i = 0;
    if Vcheck == 1
        % check for empty structure
        for i = 1:length(names_i)
            count_i = count_i + isempty(getfield(D.inv{val}.inverse,names_i{i}));
        end
        if count_i ~= length(names_i)
            disp('Inverse structure is not empty - no copy has been made for those fields');
        end
    end
    if Vcheck == 0 | count_i == length(names_i)
        % copy fields according to options
        if (Cflags(4) == 1) | (Cflags(4) == 2)
            D.inv{val}.inverse.activity     = D.inv{val-1}.inverse.activity;
            D.inv{val}.inverse.contrast     = D.inv{val-1}.inverse.contrast;
            D.inv{val}.inverse.woi          = D.inv{val-1}.inverse.woi;
        end
        if Cflags(4) == 2
            D.inv{val}.inverse.dim          = D.inv{val-1}.inverse.dim;
            D.inv{val}.inverse.priors       = D.inv{val-1}.inverse.priors;
        end
        disp('Inverse fields copied');
    end
end

% headmodel to copy for ECD
if any(Cflags) & strcmp(D.inv{val}.method,'ECD')
    if length(D.inv{val-1}.model)
        D.inv{val}.model = D.inv{val-1}.model;
    end
end

save(fullfile(D.path,D.fname),'D');
