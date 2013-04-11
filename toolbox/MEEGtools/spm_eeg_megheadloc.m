function D = spm_eeg_megheadloc(S)
% Use head localization of CTF to select/reject trials based on head
% position and (optionally) correct the sensor coordinates to correspond to
% the selected trials. The function can be used on a single dataset as well
% as several datasets together. Most of the functionality requires the
% original CTF header (read with CTF toolbox) to be present (set
% S.saveorigheader = 1 at conversion).
%
% FORMAT D = spm_eeg_megheadloc(S)
%
% S         - struct (optional)
% (optional) fields of S:
% S.D - meeg object, filename or a list of filenames of SPM EEG files
% S.rejectbetween - reject trials based on their difference from other
%                   trials (1 - yes, 0- no).
% S.threshold -     distance threshold for rejection (in meters), default
%                   0.01 (1 cm)
% S.rejectwithin -  reject trials based on excessive movement within trial
% S.trialthresh -   distance threshold for rejection (in meters), default
%                   0.005 (0.5 cm)
% S.losttrack   -   how to handle segments where the system lost track of
%                   one of the coils.
%                   'reject' - reject the trial
%                   'preserve' - try to preserve the trials. The exact
%                   behavior depends on 'rejectbetween' and 'rejectwithin'
%                   settings
%
% S.correctsens -   calculate corrected sensor representation and put it in
%                   the dataset.
% S.trialind    -   only look at a subset of trials specified. Can be used
%                   to work trial-by trial with a single file.
% S.save        -   save the header files (otherwise just return the headers).
% S.toplot      -   plot feedback information (default 1, yes).
%
% Output:
% D         - MEEG data struct or cell array of MEEG objects with the
%             rejected trials set to bad and sensors corrected (if
%             requested).
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2008 Institute of Neurology, UCL

% Vladimir Litvak, Robert Oostenveld
% $Id: spm_eeg_megheadloc.m 5396 2013-04-11 13:38:24Z vladimir $


[Finter,Fgraph,CmdLine] = spm('FnUIsetup','MEG head locations',0);

if nargin == 0
    S = [];
end

try
    D = S.D;
catch
    D = spm_select([1 inf], '\.mat$', 'Select M/EEG mat file(s)');
    S.D = D;
end

iswork = 1;
if ~isa(D, 'meeg')
    try
        for i = 1:size(D, 1)
            F{i} = spm_eeg_load(deblank(D(i, :)));
            if ~ft_senstype(chanlabels(F{i}), 'ctf')
                warning('Head localization is not supported for this data')
                iswork =0;
            end
        end
        D = F;
    catch
        error('Trouble reading files');
    end
else
    D = {D};
end
if ~iswork
    return;
end

if ~isfield(S, 'rejectbetween')
    S.rejectbetween = spm_input('Reject for intertrial difference?','+1','yes|no',[1 0], 1);
end

if ~isfield(S, 'threshold')  &&  S.rejectbetween
    S.threshold = spm_input('Threshold (cm):', '+1', 'r', '1', 1)/100;
end

if ~isfield(S, 'rejectwithin')
    S.rejectwithin = spm_input('Reject for intrial differences?','+1','yes|no',[1 0], 1);
end


if ~isfield(S, 'trialthresh') && S.rejectwithin
    S.trialthresh = spm_input('Threshold (cm):', '+1', 'r', '0.5', 1)/100;
end

if ~isfield(S, 'correctsens')
    S.correctsens = spm_input('Correct sensors?', '+1','yes|no',[1 0], 1);
end

if ~isfield(S, 'losttrack')
    S.losttrack = spm_input('How to handle losses of tracking?','+1','reject|preserve', strvcat('reject', 'preserve'));
end

if ~isfield(S, 'save')
    S.save = spm_input('Save file(s)?', '+1','yes|no',[1 0], 0);
end

if ~isfield(S, 'toplot'),         S.toplot = 1;                     end


%read HLC-channels
%HLC0011 HLC0012 HLC0013 x, y, z coordinates of nasion-coil in m.
%HLC0021 HLC0022 HLC0023 x, y, z coordinates of lpa-coil in m.
%HLC0031 HLC0032 HLC0033 x, y, z coordinates of rpa-coil in m.
hlc_chan_label = {'HLC0011' 'HLC0012' 'HLC0013'...
    'HLC0021' 'HLC0022' 'HLC0023'...
    'HLC0031' 'HLC0032' 'HLC0033'};

Ntrls=0;
dat  = [];
trlind = [];
fileind=[];

if S.toplot
    pntfig = spm_figure('GetWin','Graphics'); figure(pntfig); clf
    colors = {'b', 'g' , 'c', 'm', 'y', 'k'};
end

for f=1:numel(D)
    hlc_chan_ind = spm_match_str(chanlabels(D{f}), hlc_chan_label);
    
    Ntrl = D{f}.ntrials;
    Ntrls=Ntrls+Ntrl;
    
    if length(hlc_chan_ind) == 9
        
        if isfield(S, 'trialind') && ~isempty(S.trialind)
            trlsel = S.trialind;
        else
            trlsel = 1:Ntrl;
        end
        
        for k = trlsel
            tmpdat  = D{f}(hlc_chan_ind, :, k);
            
            header_fid = 0.01*D{f}.origheader.hc.dewar';
            cont_fid  = permute(reshape(tmpdat', [], 3, 3), [1 3 2]);
            
            dist_dev = [
                (squeeze(sqrt(sum((cont_fid(:, 1, :) - cont_fid(:, 2, :)).^2, 3))) - norm(header_fid(1,:) - header_fid(2,:)))';...
                (squeeze(sqrt(sum((cont_fid(:, 2, :) - cont_fid(:, 3, :)).^2, 3))) - norm(header_fid(2,:) - header_fid(3,:)))';...
                (squeeze(sqrt(sum((cont_fid(:, 3, :) - cont_fid(:, 1, :)).^2, 3))) - norm(header_fid(3,:) - header_fid(1,:)))'];
            
            tracking_lost_ind = find(any(abs(dist_dev) > 0.02));
            
            if ~isempty(tracking_lost_ind)
                warning(['Tracking loss detected in file ' D{f}.fname ' trial ' num2str(k)]);
                tracking_lost = 1;
                if isequal(S.losttrack, 'preserve')
                    tmpdat(:, tracking_lost_ind) = [];
                    if ~isempty(tmpdat)
                        tracking_lost = 0;
                    end
                end
            else
                tracking_lost = 0;
            end
            
            utmpdat = unique(tmpdat', 'rows')';
            
            if S.rejectwithin && ~tracking_lost
                if size(utmpdat, 2) == 1
                    dN = 0;
                    dL = 0;
                    dR = 0;
                else
                    try
                        pdist([0;1]);
                        pdistworks = 1;
                    catch
                        pdistworks = 0;
                    end
                    
                    if pdistworks
                        dN=max(pdist(utmpdat(1:3, :)'));
                        dL=max(pdist(utmpdat(4:6, :)'));
                        dR=max(pdist(utmpdat(7:9, :)'));
                    else
                        dN=max(slowpdist(utmpdat(1:3, :)'));
                        dL=max(slowpdist(utmpdat(4:6, :)'));
                        dR=max(slowpdist(utmpdat(7:9, :)'));
                    end
                end
            end
            
            if ~tracking_lost && (~S.rejectwithin || (max([dN dL dR])<S.trialthresh))
                dat     = [dat median(tmpdat, 2)];
                trlind = [trlind k];
                fileind= [fileind f];
            elseif (tracking_lost && isequal(S.losttrack, 'preserve'))
                dat     = [dat nan(9, 1)];
                trlind = [trlind k];
                fileind= [fileind f];
            else
                D{f} = badtrials(D{f}, k, 1);
            end
        end
        
        if S.toplot
            subplot('position',[0.05 0.65 0.4 0.3]);
            pltdat = squeeze(mean(reshape(tmpdat', [], 3, 3), 1))';
            
            pltdat = [pltdat; pltdat(1, :)];
            
            if ~badtrials(D{f}, k)
                h = plot3(pltdat(:, 1), pltdat(:, 2), pltdat(:, 3), ...
                    [colors{mod(f, length(colors))+1}]);
            else
                h = plot3(pltdat(:, 1), pltdat(:, 2), pltdat(:, 3), 'r');
            end
            
            hold on
            
            if isfield(D{f}, 'origheader')
                pltdat = header_fid;
            end
        end
    else
        warning(['The 9 headloc channels were not found in dataset ' D{f}.fname '. Using a single location']);
        
        % This allows for handling also files without continuous head
        % localization in a consistent way. The single location is
        % replicated according to the number of trials in the file
        
        if isfield(D{f}, 'origheader')
            utmpdat = 0.01*D{f}.origheader.hc.dewar(:);
        else
            error('Original header is required for this functionality.');
        end
        
        if S.toplot
            subplot('position',[0.05 0.65 0.4 0.3]);
            pltdat = 0.01*D{f}.origheader.hc.dewar';
        end
        
        dat     = [dat repmat(utmpdat, 1, Ntrl)];
        trlind = [trlind 1:Ntrl];
        fileind= [fileind f*ones(1, Ntrl)];
    end
    
    if S.toplot
        h = plot3(pltdat(:, 1), pltdat(:, 2), pltdat(:, 3), ...
            ['o' colors{mod(f, length(colors))+1}], 'MarkerSize', 10);
        hold on
    end
end

if S.toplot
    axis auto
    grid on
    axis equal
    axis vis3d
end

disp(['Accepted ' num2str(length(trlind)) '/' num2str(Ntrls) ' trials.']);

captured = [];

if  length(trlind)>1 && ~all(all(isnan(dat)))
    % If there was loss of tracking for just some of the trials, reject
    % them and continue with the rest.
    nanind = find(any(isnan(dat)));
    if ~isempty(nanind)
        if S.rejectbetween
            for i = length(nanind)
                D{fileind(nanind(i))} = badtrials(D{fileind(nanind(i))}, trlind(nanind(i)), 1);
            end
        end
        dat(:, nanind) = [];
        trlind(nanind) = [];
        fileind(nanind)= [];
    end
    
    if S.rejectbetween
        %%
        % Here the idea is to put a 'sphere' or 'hypercylinder' in the space of trial location whose
        % radius is 'threshold' and which captures as many trials as possible. For this we first look for the point around which the
        % density of trials is maximal. We put the cylinder there and then try to
        % optimize its position further to include more trials if possible.
        
        % The density is compute in PCA space of at most 3 dimensions
        [coeff, score, eigv] = princomp(dat');
        %%
        boundL=min(score);
        boundU=max(score);
        dim = max(find(abs(boundU-boundL)>S.threshold));
        dim=min(dim,3);
        %%
        if ~isempty(dim)
            disp(['First ' num2str(dim) ' PCs explain ' num2str(100*sum(eigv(1:dim))/sum(eigv)) '% of the variance.']);
            %%
            hcubesize=[];
            for d=1:dim
                gridres=(boundU(d)-boundL(d))./(2.^nextpow2((boundU(d)-boundL(d))/S.threshold));
                edges{d} = boundL(d):gridres:boundU(d);
                edges{d}(1)=edges{d}(1)-eps;
                edges{d}(end)=edges{d}(end)+eps;
                hcubesize=[hcubesize (length(edges{d})-1)];
            end
            
            hcube=squeeze(zeros([1 hcubesize]));
            trlpos=zeros(1, size(score,1));
            for i=1:size(score,1)
                coord=[];
                for d=1:dim
                    coord=[coord find(histc(score(i,d), edges{d}))];
                end
                
                if dim>1
                    coord=num2cell(coord);
                    trlpos(i)=sub2ind(hcubesize, coord{:});
                else
                    trlpos(i)=coord;
                end
                
                hcube(trlpos(i))=hcube(trlpos(i))+1;
            end
            %%
            [junk maxind]=max(hcube(:));
            %%
            center = mean(dat(:, find(trlpos==maxind)), 2);
        else
            disp('All trials within threshold borders');
            center = mean(dat, 2);
        end
        %%
        
        % Here the cylinder location is further optimized
        options = optimset('Display', 'iter', 'TolFun', 1, 'TolX', S.threshold/10);
        center = fminsearch(@(center) trials_captured(center, dat, S.threshold), center, options);
        %%
        % This generates the final list of trials captured in the cylinder
        [center, captured]= trials_captured(center, dat, S.threshold);
        %%
        if S.toplot
            figure(pntfig);
            subplot('position',[0.55 0.65 0.4 0.3]);
            
            plot3(score(captured,1), score(captured,2), score(captured,3),'r.');
            hold on
            plot3(score(~captured,1), score(~captured,2), score(~captured,3),'.');
            
            axis equal
        end
        %%
        ufileind=unique(fileind(captured));
        %%
        
        for f = 1:numel(D)
            origreject = badtrials(D{f});
            D{f} = badtrials(D{f}, ':', 1);
            if ismember(f, ufileind) && any(captured & (fileind == f))
                D{f} = badtrials(D{f}, trlind(captured & (fileind == f)), 0);
                if any(origreject)
                    D{f} = badtrials(D{f}, find(origreject), 1);
                end
            end
        end
    end
end

if isempty(captured)
    captured = ones(1, size(dat, 2));
end


%%
% This generates the corrected grad structure
if S.correctsens && ((length(hlc_chan_ind) == 9) || numel(D)>1) && ~isempty(trlind) && ~all(all(isnan(dat)))
    
    newcoils=mean(dat(:, captured), 2);
    
    nas = newcoils(1:3);
    lpa = newcoils(4:6);
    rpa = newcoils(7:9);
    
    %compute transformation matrix from dewar to head coordinates
    M = spm_eeg_inv_headcoordinates(nas, lpa, rpa);
    
    dewar = 0.01*D{1}.origheader.hc.dewar;
    
    M1 = spm_eeg_inv_headcoordinates(dewar(:,1)', dewar(:,2)', dewar(:,3)');
    
    grad = sensors(D{1}, 'MEG');
    
    newgrad = ft_transform_sens(M*inv(M1), grad);
    
    if S.toplot
        figure(pntfig);
        subplot('position',[0.05 0.05 0.9 0.5]);
        
        
        cfg = [];
        cfg.style = '3d';
        cfg.rotate = 0;
        cfg.grad = grad;
        
        lay = ft_prepare_layout(cfg);
        
        cfg.grad = newgrad;
        
        newlay = ft_prepare_layout(cfg);
        
        plot3(newlay.pos(:,1), newlay.pos(:,2), newlay.pos(:,3), '.r', 'MarkerSize', 5);
        
        hold on
        
        plot3(lay.pos(:,1), lay.pos(:,2), lay.pos(:,3), '.k');
        
        axis equal off
    end
    
    newfid = ft_transform_headshape(M*inv(M1), fiducials(D{1}));
    
    for f = 1:numel(D)
        D{f} = sensors(D{f}, 'MEG', newgrad);
        D{f} = fiducials(D{f}, newfid);
    end
end

for f = 1:numel(D)
    D{f} = history(D{f}, 'spm_eeg_megheadloc', S);
    if S.save
        save(D{f});
    end
end

if numel(D) == 1
    D = D{1};
end

end

function [obj, captured] = trials_captured(center, dat, threshold)

nsmp    = size(dat,2);
dist = max(sqrt([1 1 1 0 0 0 0 0 0; 0 0 0 1 1 1 0 0 0; 0 0 0 0 0 0 1 1 1]*(dat - repmat(center(:), 1, nsmp)).^2));
obj = -sum(dist<threshold);

if nargout>1
    captured=dist<threshold;
end

end

function Y = slowpdist(X)
% Thanks to Guillaume
N = size(X,1);
Y = zeros(1,N*(N-1)/2);
k = 1;
for i=1:N-1
    Y(k:(k+N-i-1)) = sqrt(sum((repmat(X(i,:),N-i,1) - X((i+1):N,:)).^2,2));
    k = k + N - i;
end
end
