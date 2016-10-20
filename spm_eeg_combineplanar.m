function D = spm_eeg_combineplanar(S)
% Combines data from MEGPLANAR sensors
% FORMAT D = spm_eeg_combineplanar(S)
%
% S        - optional input struct
%  fields of S:
%   D        - MEEG object or filename
%   mode     -
%              'append'     - add combined channels to the origal channels
%              'replace'    - replace MEGPLANAR with combined (default)
%              'replacemeg' - replace all MEG channels with combined but
%                             keep non-MEG
%              'keep'       - only write out the combined channels
%
%   prefix   - prefix for the output file (default - 'P')
%
%
% Output:
% D        - MEEG object (also written on disk)
%
%__________________________________________________________________________
% Copyright (C) 2008-2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_combineplanar.m 6904 2016-10-20 12:04:59Z vladimir $

SVNrev = '$Rev: 6904 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','MEG Combine planar'); spm('Pointer','Watch');

if ~isfield(S, 'prefix'),       S.prefix   = 'P';           end
if ~isfield(S, 'mode'),         S.mode     = 'replace';     end

D = spm_eeg_load(S.D);

isTF = strncmpi(D.transformtype,'TF',2);

chanset = spm_eeg_planarchannelset(D.chanlabels);

chanind = [];
labelnew = {};
for i = 1:size(chanset, 1)
    cind = D.indchannel(chanset(i, 1:2));
    if length(cind) == 2
        chanind  = [chanind cind];
        labelnew = [labelnew; chanset(i, 4)];
    end
end

megind    = D.indchantype('MEG');
planarind = D.indchantype('MEGPLANAR');

% add a row of zeros;
chanind(3, end) = 0;
copyind = [];
switch S.mode
    case 'append'
        copyind = 1:D.nchannels;
        copyind = [copyind; copyind];
        chanind(3, :) = (D.nchannels+1):(D.nchannels+size(chanind, 2));
    case {'replace', 'replacemeg'}
        ind = 1;
        for i = 1:D.nchannels
            if any(i == planarind)
                k = find((chanind(1, :) == i) | (chanind(2, :) == i));
                if ~isempty(k) && (chanind(3, k) == 0)
                    chanind(3, k) = ind;
                    ind = ind + 1;
                end
            elseif isequal(S.mode, 'replace') || ~any(i == megind)
                copyind = [copyind [i ind]'];
                ind = ind + 1;
            end
        end
    case 'keep'
        copyind       = [];
        chanind(3, :) = 1:size(chanind, 2);
end

if isempty(copyind)
    Nchannels = max(chanind(3, :));
else
    Nchannels = max([chanind(3, :), copyind(2, :)]);
end

%-Generate new MEEG object with new files
%--------------------------------------------------------------------------
if isTF
    Dnew = clone(D, [S.prefix fname(D)], [Nchannels D.nfrequencies D.nsamples D.ntrials], 1);
else
    Dnew = clone(D, [S.prefix fname(D)], [Nchannels D.nsamples D.ntrials], 1);
end

if strcmp(D.type, 'continuous')
    blksz  = D.fsample;
    blknum = floor(D.nsamples/blksz);
    
    %----------------------------------------------------------------------
    spm_progress_bar('Init', blknum, 'Data blocks processed');
    if blknum > 100, Ibar = floor(linspace(1, blknum,100));
    else Ibar = 1:blknum; end
    
    if isTF
        i = 1;
        while 1
            if i<=blknum
                Iblock  = ((i-1)*blksz + 1):(i*blksz);
            elseif (i == (blknum + 1)) && blknum*blksz < D.nsamples
                Iblock  = (blknum*blksz + 1):D.nsamples;
            else
                break;
            end
            
            
            planar1 = D(chanind(1, :), :, Iblock);
            planar2 = D(chanind(2, :), :, Iblock);
                       
            Dnew(chanind(3, :), :, Iblock) = planar1 + planar2;
            
            
            if ~isempty(copyind)
                Dnew(copyind(2, :), :, Iblock) =  D(copyind(1, :), :, Iblock);
            end
            
            if any(Ibar == i), spm_progress_bar('Set', i); end
            
            i = i+1;
        end
    else
        i = 1;
        while 1
            if i<=blknum
                Iblock  = ((i-1)*blksz + 1):(i*blksz);
            elseif (i == (blknum + 1)) && blknum*blksz < D.nsamples
                Iblock  = (blknum*blksz + 1):D.nsamples;               
            else
                break;
            end
            
            
            planar1 = D(chanind(1, :),  Iblock);
            planar2 = D(chanind(2, :),  Iblock);
            
            Dnew(chanind(3, :), Iblock) = sqrt(planar1.^2 + planar2.^2);
            
            if ~isempty(copyind)
                Dnew(copyind(2, :), Iblock) =  D(copyind(1, :), Iblock);
            end
            
            if any(Ibar == i), spm_progress_bar('Set', i); end
            
            i = i+1;
        end
    end
else  %-------- Epoched data --------------------------------------------------
    spm_progress_bar('Init', D.ntrials, 'Trials processed');
    if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
    else Ibar = 1:D.ntrials; end
    
    if isTF
        for i=1:D.ntrials
            planar1 = D(chanind(1, :), :, :, i);
            planar2 = D(chanind(2, :), :, :, i);
            
            Dnew(chanind(3, :),:, :, i) = planar1 + planar2;
            
            if ~isempty(copyind)
                Dnew(copyind(2, :), :, :, i) =  D(copyind(1, :), :, :, i);
            end
            
            if any(Ibar == i), spm_progress_bar('Set', i); end
        end
    else
        for i=1:D.ntrials
            planar1 = D(chanind(1, :), :, i);
            planar2 = D(chanind(2, :), :, i);
            
            Dnew(chanind(3, :), :, i) = sqrt(planar1.^2 + planar2.^2);
            
            if ~isempty(copyind)
                Dnew(copyind(2, :), :, i) =  D(copyind(1, :), :, i);
            end
            
            if any(Ibar == i), spm_progress_bar('Set', i); end
        end
    end
end

spm_progress_bar('Clear');

Dnew = chanlabels(Dnew, chanind(3, :), labelnew);
Dnew = badchannels(Dnew, chanind(3, :), badchannels(D, chanind(1, :)) | badchannels(D, chanind(2, :)));
Dnew = chantype(Dnew, chanind(3, :), 'MEGCOMB');
Dnew = units(Dnew, chanind(3, :), units(D, chanind(1, :)));
Dnew = coor2D(Dnew, chanind(3, :), 0.5*(coor2D(D, chanind(1,:))+coor2D(D, chanind(1,:))));

if ~isempty(copyind)
    Dnew = chanlabels(Dnew, copyind(2, :), D.chanlabels(copyind(1, :)));
    Dnew = badchannels(Dnew, copyind(2, :), badchannels(D, copyind(1, :)));
    Dnew = chantype(Dnew, copyind(2, :), D.chantype(copyind(1, :)));
    Dnew = units(Dnew, copyind(2, :), D.units(copyind(1, :)));
    Dnew = coor2D(Dnew, copyind(2, :),coor2D(D, copyind(1, :)));
end

%-Save the new M/EEG dataset
%--------------------------------------------------------------------------
Dnew = Dnew.history(mfilename, S);
save(Dnew);

D = Dnew;

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','MEG Combine planar: done'); spm('Pointer','Arrow');
