function Dnew = spm_eeg_remove_spikes(S)
% Use Will Pennys robust GLM code to remove 'spikes' from continuous data.
% Such spikes occur in EEG data recorded with the CTF MEG system at FIL
% due to some obscure electrical problem.
%
% FORMAT Dnew = spm_eeg_remove_spikes(S)
%
% S         - struct (optional)
% (optional) fields of S:
% S.D      - meeg object or filename
% S.logbf  - clean a block if log bayes factor in favour of spike model is
%            bigger than this (default - 3)
% S.hpf    - high-pass frequency above which to look for spikes (default 40 Hz)
% S.fast   - option to speed up the function by only using GLM if there is
%            threshold crossing ('yes', or check all the data with GLM - 'no')
% S.fasthresh - threshold for the fast option (in STD) - default 4
% S.trialbased - use trials in the data as they are ('yes') or break them
%              into sub-blocks ('no' - default)
% S.channels - channels to clean up (default 'gui' - brings up a GUI for
%             channel choice.
%
% Output:
% Dnew  - MEEG object with data cleaned of spikes.
%
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2008 Institute of Neurology, UCL

% Vladimir Litvak, Will Penny
% $Id: spm_eeg_remove_spikes.m 5640 2013-09-18 12:02:29Z vladimir $

if nargin == 0
    S = [];
end

% clean a block if log bayes factor in favour of spike model is bigger than this
if ~isfield(S, 'logbf'),             S.logbf = 3;                            end
% Have signals up to this freq in design matrix
if ~isfield(S, 'hpf'),               S.hpf = 40;                             end
% Use 1s blocks
if ~isfield(S, 'fast'),              S.fast='yes';                           end
if ~isfield(S, 'fasthresh'),         S.fasthresh=4;                          end
if ~isfield(S, 'trialbased'),        S.trialbased='no';                      end
if ~isfield(S, 'channels'),          S.channels='gui';                       end

try
    D = S.D;
catch
    [D,sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

if ~isfield(S, 'blocksize'),         S.blocksize = D.fsample;                end

spm('Pointer', 'Watch');drawnow;

if ~isfield(S, 'channels')
    chansel = D.indchantype({'MEEG', 'LFP'});
    S.channels = D.chanlabels(chansel);
elseif ischar(S.channels) && strcmp(S.channels, 'gui')
    chansel = listdlg('ListString', D.chanlabels, 'SelectionMode', 'multiple' ,'Name', 'Select channels' , 'ListSize', [400 300]);
    S.channels = D.chanlabels(chansel);
else
    chansel     = spm_match_str(D.chanlabels, S.channels);
end

nchan = length(chansel);    % number of channels

if strcmpi(S.trialbased, 'no') || (D.ntrials == 1)
    Nblocks=floor(D.nsamples./S.blocksize);
    block_size=S.blocksize;
    trialbased=0;
else
    Nblocks=D.ntrials;
    block_size=D.nsamples;
    trialbased=1;
end

% Make design matrix
K(1).row=[1:block_size];
K(1).RT=1/D.fsample;
K(1).HParam=1/S.hpf;
K = spm_filter(K);
X=K.X0;
X=[X ,ones(size(X,1),1)];


if ~trialbased && (D.ntrials > 1)
    rptnum = D.ntrials;
else
    rptnum = 1;
end

% generate new meeg object with new filenames
Dnew = clone(D, ['S' fnamedat(D)], [D.nchannels D.nsamples D.ntrials]);

for i = 1:D.nchannels
    Dnew(i, :, :) = D(i, :, :);
end

for r = 1:rptnum
    for ch=1:nchan
        for n=1:Nblocks,
            if trialbased
                lfp=squeeze(D(chansel(ch), :, n));
            else
                si=(n-1)*block_size+1;
                ei=si+block_size-1;
                lfp=squeeze(D(chansel(ch), si:ei, r));
            end
            m_lfp=mean(lfp);
            s_lfp=std(lfp);
            if s_lfp==0
                continue;
            end
            
            lfp=(lfp-m_lfp)/s_lfp;
            
            if ~(strcmpi(S.fast, 'yes') && max(abs(lfp))<S.fasthresh)
                glm=spm_rglm(lfp,X,1);
                [rglm,lfp_clean]=spm_rglm(lfp,X,2);
                if rglm.fm-S.logbf > glm.fm
                    if trialbased
                        Dnew(chansel(ch), :, n)=s_lfp*lfp_clean'+m_lfp;
                    else
                        Dnew(chansel(ch), si:ei, r) = s_lfp*lfp_clean'+m_lfp;
                    end
                    descr='spiky';
                else
                    descr='clean';
                end
                disp(sprintf('Channel %d out of %d, trial %d out of %d, block %d out of %d %s',ch,nchan, r, rptnum, n,Nblocks,descr));
            end
        end
    end
end

Dnew = history(Dnew, 'spm_eeg_remove_spikes', S);

save(Dnew);

spm('Pointer', 'Arrow');drawnow;
