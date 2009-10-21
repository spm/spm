function H = spm_eeg_history(S)
% Generate a MATLAB script from the history of an M/EEG SPM data file
% FORMAT H = spm_eeg_history(S)
%
% S  - filename or input struct (optional)
% (optional) fields of S:
% history         - history of M/EEG object (D.history)
% sname           - filename of the to be generated MATLAB script
%
% H               - cell array summary of history for review purposes
%__________________________________________________________________________
%
% In SPM for M/EEG, each preprocessing step enters its call and input
% arguments into an internal history. The sequence of function calls that
% led to a given file can be read by the history method (i.e. call
% 'D.history'). From this history this function generates a script (m-file)
% which can be run without user interaction and will repeat, if run, the
% exact sequence on the preprocessing steps stored in the history. Of
% course, the generated script can also be used as a template for a
% slightly different analysis or for different subjects.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_history.m 3497 2009-10-21 21:54:28Z vladimir $

try
    h = S.history;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, H = {}; return; end
    D = spm_eeg_load(D);
    h = D.history;
end

if nargout
    H = hist2cell(h);
else
    try
        S.sname; % not fname, as it means something else for MEEG object
    catch
        [filename, pathname] = uiputfile('*.m', ...
            'Select the to be generated script file');
        if isequal(filename,0) || isequal(pathname,0)
            return;
        end
        S.sname = fullfile(pathname, filename);
    end

    hist2script(h,S.sname);
end

    
%==========================================================================
% function hist2script
%==========================================================================
function hist2script(h,fname)

histlist = convert2humanreadable(h);
[selection, ok]= listdlg('ListString', histlist, 'SelectionMode', 'multiple',...
    'InitialValue', 1:numel(histlist) ,'Name', 'Select history entries', ...
    'ListSize', [400 300]);
if ok, h = h(selection); else return; end

Nh = length(h);
fp = fopen(fname, 'wt');
if fp == -1
    error('File %s cannot be opened for writing.', fname);
end

fprintf(fp, '%s\n\n', 'spm(''defaults'', ''eeg'');');

for i = 1:Nh

    fprintf(fp, '%s\n', 'S = [];');
    
    s = gencode(h(i).args(1), 'S');

    for j = 1:length(s)
        fprintf(fp, '%s\n', s{j});
    end
    
    fprintf(fp, '%s\n\n\n', ['D = ' h(i).fun '(S);']);
    
end

fclose(fp);


%==========================================================================
% function hist2cell
%==========================================================================
function H = hist2cell(h)

nf = length(h);
H = cell(nf,4);
for i=1:nf
    H{i,1} = convert2humanreadable(h(i));
    H{i,2} = h(i).fun;
    args = h(i).args;
    try,args = args{1};end
    switch h(i).fun
        case 'spm_eeg_convert'
            H{i,3} = args.dataset;
            if i<nf
                path = fileparts(H{i,3});
                H{i,4} = fullfile(path,[args.outfile '.mat']);
            else
                H{i,4} = '[this file]';
            end
        case 'spm_eeg_prep'
            Df = args.D;
            try, Df = args.D.fname; end
            H{i,2}  = [H{i,2},' (',args.task,')'];
            pth     = fileparts(Df);
            if isempty(pth)
                try
                    pth = fileparts(H{i-1,4});
                    Df  = fullfile(pth, Df);
                end
            end
            H{i,3} = Df;
            if i<nf
                H{i,4} = Df;
            else
                H{i,4} = '[this file]';
            end
        otherwise
            Df = args.D;
            try,Df = args.D.fname;end
            H{i,3} = Df;
            if i<nf
                try
                    args2 = h(i+1).args;
                    try,args2 = args2{1};end
                    Df2 = args2.D;
                    try,Df2 = args2.D.fname;end
                    H{i,4} = Df2;
                catch
                    H{i,4} = '?';
                end
            else
                H{i,4} = '[this file]';
            end
    end
end


%==========================================================================
% function convert2humanreadable
%==========================================================================
function hh = convert2humanreadable(h)

hh = cell(numel(h),1);
for i=1:numel(h)
    switch h(i).fun
        case 'spm_eeg_convert'
            hh{i} = 'Convert';
        case 'spm_eeg_epochs'
            hh{i} = 'Epoch';
        case 'spm_eeg_filter'
            if length(h(i).args.filter.PHz) == 2
                hh{i} = [upper(h(i).args.filter.band(1)) h(i).args.filter.band(2:end)...
                    ' filter ' num2str(h(i).args.filter.PHz(:)', '%g %g') ' Hz'];
            else
                hh{i} = [upper(h(i).args.filter.band(1)) h(i).args.filter.band(2:end)...
                    ' filter ' num2str(h(i).args.filter.PHz, '%g') ' Hz'];
            end
        case 'spm_eeg_downsample'
            hh{i}  = ['Downsample to ' num2str(h(i).args.fsample_new) ' Hz'];
        case 'spm_eeg_bc'
            hh{i} = ['Baseline correction ' mat2str(h(i).args.time(:)') ' ms'];
        case 'spm_eeg_copy'   
            hh{i} = 'Copy dataset';
        case 'spm_eeg_montage'
            hh{i} = 'Change montage';
        case 'spm_eeg_artefact'
            hh{i} = 'Detect artefacts';
        case 'spm_eeg_average'
            hh{i} = 'Average';
        case 'spm_eeg_average_TF'
            hh{i} = 'Average time-frequency';
        case 'spm_eeg_grandmean'
            hh{i} = 'Grand mean';
        case 'spm_eeg_merge'
            hh{i} = 'Merge';
        case 'spm_eeg_tf'
            hh{i} = 'Compute time-frequency';
        case 'spm_eeg_weight_epochs'
            hh{i} = 'Compute contrast';
        case 'spm_eeg_sort_conditions'
            hh{i} = 'Sort conditions';
        case 'spm_eeg_prep'
            switch h(i).args.task
                case 'settype'
                    hh{i} = 'Set channel type';
                case {'loadtemplate', 'setcoor2d', 'project3D'}
                    hh{i} = 'Set 2D coordinates';
                case 'loadeegsens'
                    hh{i} = 'Load EEG sensor locations';
                case 'defaulteegsens'
                    hh{i} = 'Set EEG sensor locations to default';
                case 'sens2chan'
                    hh{i} = 'Specify initial montage';
                case 'headshape'
                    hh{i} = 'Load fiducials/headshape';
                case 'coregister'
                    hh{i} = 'Coregister';
                otherwise
                    hh{i} = ['Prepare: ' h(i).args.task];
            end
        otherwise
            hh{i} = h(i).fun;
    end
end
if numel(h) == 1, hh = hh{1}; end
