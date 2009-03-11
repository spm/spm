function [D, S] = spm_eeg_convert2images(S)
% User interface for conversion of M/EEG-files to SPM image file format
% FORMAT [D, S] = spm_eeg_convert2images(S)
%
% S                   - input structure (optional)
% (optional) fields of S:
%   S.D               - MEEG object or filename of M/EEG mat-file with
%                       epoched data
%   S.images with entries (all optional):
%     fmt             - string that determines type of input file. Currently,
%                       it can be 'channels' or 'frequency'
%     elecs           - electrodes of interest (as vector of indices)
%     region_no       - region number
%     freqs           - frequency window of interest (2-vector) [Hz]
%   S.n               - dimension of output images in voxels (scalar because
%                       output will be square image)
%   S.interpolate_bad - flag (0/1) that indicates whether values for
%                       bad channels should be interpolated (1) or left
%                       out (0).
% output:
% S         - can be used to construct script (as in the history-function)
%__________________________________________________________________________
%
% spm_eeg_convert2images is a user interface to convert M/EEG-files in SPM
% format to SPM's image format, using an interpolation on a 2D-plane.
% This function assembles some necessary information before branching to
% the format-specific conversion routines.
%
% Output: The converted data are written to files. The header structs, but
% not the data, are returned in D as a cell vector of structs, and the
% struct S is returned to allow for saving the history of function calls.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% James Kilner, Stefan Kiebel
% $Id: spm_eeg_convert2images.m 2857 2009-03-11 13:21:04Z guillaume $

SVNrev = '$Rev: 2857 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG conversion setup'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat$', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

if strcmp(D.type, 'continuous')
    error('Data are continuous. Try epoched data.');
end

%-Time-Frequency data
%==========================================================================
if strncmpi(D.transformtype, 'TF',2);

    %-Average over channels or frequencies?
    %----------------------------------------------------------------------
    try
        images.fmt = S.images.fmt;
    catch
        spm_input('average over ...', 1, 'd')
        Ctype        = {'channels', 'frequency'};
        str          = 'Average over which dimension';
        Sel          = spm_input(str, 2, 'm', Ctype);
        images.fmt   = Ctype{Sel};
        S.images.fmt = images.fmt;
    end

    switch images.fmt
        %-Average over channels
        %------------------------------------------------------------------
        case {'channels'}

            %-Select channels
            %--------------------------------------------------------------
            try
                images.electrodes_of_interest = S.images.elecs;
            catch
                str  = 'electrodes[s]';
                Ypos = '+1';
                while true
                    [images.electrodes_of_interest, Ypos] = ...
                        spm_input(str, Ypos, 'r', [], [1 Inf]);
                    if any(ismember(images.electrodes_of_interest, [1:D.nchannels]))
                        break;
                    end
                end
                S.images.elecs = images.electrodes_of_interest;
            end

            %-Attribute a region number to that channel
            %--------------------------------------------------------------
            try
                images.Nregion = S.images.region_no;
            catch
                str = 'region number';
                images.Nregion = spm_input(str, '+1', 'r', [], [1 Inf]);
                S.images.region_no = images.Nregion;
            end

            %-Convert to NIfTI images
            %--------------------------------------------------------------
            cl  = D.condlist;
            for i = 1 : D.nconditions
                Itrials = pickconditions(D, cl{i}, 1)';

                dname = sprintf('%dROI_TF_trialtype_%s', images.Nregion, cl{i});
                [sts, msg] = mkdir(D.path, dname);
                if ~sts, error(msg); end
                P = fullfile(D.path, dname);

                for l = Itrials(:)'

                    if strcmp(D.type, 'single')
                        % single trial data
                        fname = sprintf('trial%04d.img', l);
                    else
                        % evoked data
                        fname = 'average.img';
                    end
                    fname = fullfile(P,fname);

                    N     = nifti;
                    dat   = file_array(fname, [D.nfrequencies D.nsamples], 'FLOAT64-LE');
                    N.dat = dat;
                    N.mat = [...
                        1  0               0  min(D.frequencies);...
                        0  1000/D.fsample  0  time(D, 1, 'ms');...
                        0  0               1  0;...
                        0  0               0  1];
                    N.mat_intent = 'Aligned';
                    create(N);

                    N.dat(:, :) = spm_cond_units(squeeze(mean(D(images.electrodes_of_interest, :, :, l), 1)));

                end
            end

        %-Average over frequency
        %------------------------------------------------------------------
        case {'frequency'}

            %-Select frequency window
            %--------------------------------------------------------------
            try
                images.Frequency_window = S.images.freqs;
                inds = find(D.frequencies >= images.Frequency_window(1) & ...
                    D.frequencies <= images.Frequency_window(2));
            catch
                str = 'Frequency window';
                Ypos = '+1';
                while true
                    [images.Frequency_window, Ypos] = spm_input(str, Ypos, 'r', [], 2);
                    inds = find(D.frequencies >= images.Frequency_window(1) & ...
                        D.frequencies <= images.Frequency_window(2));
                    if ~isempty(inds), break; end
                    str = 'No data in range';
                end
                S.images.freqs = images.Frequency_window;
            end

            %-Generate new dataset with averaged data over frequency window
            %--------------------------------------------------------------
            fnamedat = ['F' num2str(images.Frequency_window(1)) '_' ...
                num2str(images.Frequency_window(2)) '_' D.fnamedat];
            Dnew = clone(D, fnamedat, [D.nchannels D.nsamples D.ntrials]);

            Dnew(1:Dnew.nchannels, 1:Dnew.nsamples, 1:Dnew.ntrials) = ...
                squeeze(mean(D(:, inds, :, :), 2));
            Dnew = transformtype(Dnew, 'time');
            save(Dnew);

            %-Convert that dataset into images
            %--------------------------------------------------------------
            S.Fname = fullfile(Dnew.path, Dnew.fname);
            S = spm_eeg_convert2scalp(S);

        %-Otherwise...
        %------------------------------------------------------------------
        otherwise
            error('Unknown dimension to average over.');
    end
else

    %-Time Epoched data
    %======================================================================
    S.Fname = fullfile(D.path, D.fname);
    S = spm_eeg_convert2scalp(S);
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG conversion: done'); spm('Pointer','Arrow');
