function [D, S, Pout] = spm_eeg_convert2images(S)
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
%     elecs           - channels of interest (as vector of indices)
%     region_no       - region number
%     freqs           - frequency window of interest (2-vector) [Hz]
%     t_win           - [t1 t2] For 'frequency' option with TF data, specify
%                        this field to only extract power in restricted
%                        time window. This allows you to avoid eg. edge
%                        effects.
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
% $Id: spm_eeg_convert2images.m 3724 2010-02-16 12:16:57Z vladimir $

SVNrev = '$Rev: 3724 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG conversion setup'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

if strcmp(D.type, 'continuous')
    error('Data are continuous. Try epoched data.');
end

%-Time-Frequency data
%==========================================================================
if strncmpi(D.transformtype, 'TF',2)
    
    %-If it's clear what to average over, assign automatically
    %----------------------------------------------------------------------
    if D.nfrequencies == 1
        S.images.fmt = 'frequency';
    end
    
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
                images.channels_of_interest = S.images.elecs;
            catch
                if D.nchannels > 1
                    meegchan = D.meegchannels;
                    [selection, ok]= listdlg('ListString', D.chanlabels(meegchan), 'SelectionMode', 'multiple' ,'Name', 'Select channels' , 'ListSize', [400 300]);
                    if ~ok
                        return;
                    end
                    
                    images.channels_of_interest = meegchan(selection);
                    
                else
                    images.channels_of_interest = 1;
                end
                S.images.elecs = images.channels_of_interest;
            end

            %-Attribute a region number to those channels
            %--------------------------------------------------------------
            try
                images.Nregion = S.images.region_no;
            catch
                if D.nchannels > 1
                    str = 'region number';
                    images.Nregion = spm_input(str, '+1', 'r', [], [1 Inf]);
                    S.images.region_no = images.Nregion;
                else
                    images.Nregion = 1;
                end
            end

            %-Convert to NIfTI images
            %--------------------------------------------------------------
            cl  = D.condlist;
            df  = diff(D.frequencies);
            if any(diff(df))
                warning('Irregular frequency spacing.');
            end
            
  %-Make output directory for each dataset
  %--------------------------------------------------------------------------
            [P, F] = fileparts(S.D);
            if isempty(P), P = pwd; end
            [sts, msg] = mkdir(P, F);
            if ~sts, error(msg); end
            P  = fullfile(P, F);

            Pout = cell(1, D.nconditions);
            for i = 1 : D.nconditions
                Itrials = pickconditions(D, cl{i}, 1)';                
                
                Pout{i} = {};
                
  %-Make subdirectory for each condition
  %--------------------------------------------------------------------------
                dname = sprintf('%dROI_TF_trialtype_%s', images.Nregion, cl{i});
                [sts, msg] = mkdir(P, dname);
                if ~sts, error(msg); end
                dname = fullfile(P, dname);

                for l = Itrials(:)'

                    if strcmp(D.type, 'single')
                        % single trial data
                        fname = sprintf('trial%04d.img', l);
                    else
                        % evoked data
                        fname = 'average.img';
                    end
                    fname = fullfile(dname,fname);
                    
                    Pout{i} = [Pout{i}, {fname}];

                    N     = nifti;
                    dat   = file_array(fname, [D.nfrequencies D.nsamples], 'FLOAT64-LE');
                    N.dat = dat;
                    N.mat = [...
                        df(1)   0               0  min(D.frequencies);...
                        0       1000/D.fsample  0  time(D, 1, 'ms');...
                        0       0               1  0;...
                        0       0               0  1];
                    N.mat(1,4) = N.mat(1,4) - N.mat(1,1);
                    N.mat(2,4) = N.mat(2,4) - N.mat(2,2);
                    N.mat_intent = 'Aligned';
                    create(N);

                    if ~isempty(strmatch('MEG', D.chantype(images.channels_of_interest))) &&...
                            isempty(strmatch('fT', D.units(images.channels_of_interest))) && ...
                            isempty(strmatch('dB', D.units(images.channels_of_interest))) && ...
                            isempty(strmatch('%', D.units(images.channels_of_interest)))
                        scale = 1e30;
                    else
                        scale = 1;
                    end
                    
                    N.dat(:, :) = scale*spm_squeeze(mean(D(images.channels_of_interest, :, :, l), 1), 1);

                end
                Pout{i} = char(Pout{i});
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

            % This is a slightly ugly fix for the problem of very small
            % power values for MEG (in T^2). 
            megchanind = strmatch('MEG', D.chantype);
            nonmegchanind = setdiff(1:D.nchannels, megchanind);
            
            %-Generate new dataset with averaged data over frequency window
            %--------------------------------------------------------------
            fnamedat = ['F' num2str(images.Frequency_window(1)) '_' ...
                num2str(images.Frequency_window(2)) '_' D.fnamedat];
            
            if isfield(S.images,'t_win')
                % Only extract time points in specified window
                tims=time(D);
                if S.images.t_win(1) < tims(1) || S.images.t_win(2) > tims(end)
                    disp('Error: Impossible specification of time extraction window');
                    disp(S.images.t_win);
                    return;
                end
                tind=find(tims > S.images.t_win(1) & tims < S.images.t_win(2));
                Nind=length(tind);
                Dnew = clone(D, fnamedat, [D.nchannels Nind D.ntrials]);
                
                if ~isempty(megchanind)
                    Dnew(megchanind, 1:Dnew.nsamples, 1:Dnew.ntrials) = ...
                        1e30*spm_squeeze(mean(D(megchanind, inds, tind(1):tind(end), :), 2), 2);
                end
                if ~isempty(nonmegchanind)
                    Dnew(nonmegchanind, 1:Dnew.nsamples, 1:Dnew.ntrials) = ...
                        spm_squeeze(mean(D(nonmegchanind, inds, tind(1):tind(end), :), 2), 2);
                end
                
                Dnew = timeonset(Dnew, tims(tind(1)));
            else
                Dnew = clone(D, fnamedat, [D.nchannels D.nsamples D.ntrials]);
                if ~isempty(megchanind)
                    Dnew(megchanind, 1:Dnew.nsamples, 1:Dnew.ntrials) = ...
                        1e30*spm_squeeze(mean(D(megchanind, inds, :, :), 2), 2);
                end
                if ~isempty(nonmegchanind)
                    Dnew(nonmegchanind, 1:Dnew.nsamples, 1:Dnew.ntrials) = ...
                        spm_squeeze(mean(D(nonmegchanind, inds, :, :), 2), 2);
                end
            end
                        
            Dnew = transformtype(Dnew, 'time');
            
            if ~isempty(megchanind)
                Dnew = units(Dnew, megchanind, 'fT^2');
            end

            save(Dnew);

            %-Convert that dataset into images
            %--------------------------------------------------------------
            S.Fname = fullfile(Dnew.path, Dnew.fname);
            [S, Pout] = spm_eeg_convert2scalp(S);

        %-Otherwise...
        %------------------------------------------------------------------
        otherwise
            error('Unknown dimension to average over.');
    end
else

    %-Time Epoched data
    %======================================================================
    S.Fname = fullfile(D.path, D.fname);
    [S, Pout] = spm_eeg_convert2scalp(S);
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG conversion: done'); spm('Pointer','Arrow');
