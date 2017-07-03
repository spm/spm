function D = spm_eeg_spatial_confounds(S)
% This function defines spatial confounds and adds them to MEEG dataset.
% FORMAT D = spm_eeg_spatial_confounds(S)
%
% S        - optional input struct
%  fields of S:
%   D        - MEEG object or filename of M/EEG mat-file with epoched data
%   mode     - method for definition of the confounds (EYES, BESA, SVD,
%              SPMEEG, CLEAR)
%
%
% _______________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_spatial_confounds.m 7128 2017-07-03 11:58:47Z vladimir $


SVNrev = '$Rev: 7128 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Define spatial confounds');

if ~isfield(S, 'mode') && isfield(S, 'method'),          S.mode  = S.method;             end
if ~isfield(S, 'threshold') && isfield(S, 'svdthresh'),  S.threshold  = S.svdthresh;    end


D = spm_eeg_load(S.D);

if ~isfield(S, 'conditions') || isempty(S.conditions),   S.conditions = D.condlist;      end 
if ~iscell(S.conditions), S.conditions = {S.conditions};                                 end

sconf = [];
switch upper(S.mode)
    case 'EYES'
        [D, ok] = check(D, 'sensfid');
        
        if ~ok
            if check(D, 'basic')
                errordlg(['The requested file is not ready for source reconstruction.'...
                    'Use prep to specify sensors and fiducials.']);
            else
                errordlg('The meeg file is corrupt or incomplete');
            end
            return
        end
        
        %% ============ Find or prepare head model
        
        if ~isfield(D, 'val')
            D.val = 1;
        end
        
        if ~isfield(D, 'inv') || ~iscell(D.inv) ||...
                ~(isfield(D.inv{D.val}, 'forward') && isfield(D.inv{D.val}, 'datareg')) ||...
                ~isa(D.inv{D.val}.mesh.tess_ctx, 'char') % detects old version of the struct
            D = spm_eeg_inv_mesh_ui(D, D.val);
            D = spm_eeg_inv_datareg_ui(D, D.val);
            D = spm_eeg_inv_forward_ui(D, D.val);
            
            save(D);
        end
        
        eyes = [-34 53 -38; 34 53 -38];
        
        
        sconf = [];
        sconf.label = D.chanlabels(D.indchantype('MEEG'));
        sconf.coeff = nan(length(sconf.label), 6);
        sconf.bad = ones(length(sconf.label), 1);
        
        [junk, modalities] = modality(D, 1, 1);
        
        for k = 1:numel(modalities)
            chanind = indchantype(D, modalities{k}, 'GOOD');
            
            if isempty(chanind)
                continue;
            end
            
            data = spm_eeg_inv_get_vol_sens(D, [], [], [], modalities{k});
            
           
            eyes = spm_eeg_inv_transform_points(inv(data.transforms.toMNI), eyes);
            
            vol  = data.(modalities{k}).vol;
            sens = data.(modalities{k}).sens;
            
            if isa(vol, 'char')
                vol = ft_read_vol(vol);
            end            
            
            [vol, sens] = ft_prepare_vol_sens(vol, sens, 'channel', D.chanlabels(chanind));
            
            
            if strncmp(modalities{k}, 'MEG', 3)
                reducerank = 2;
            else
                reducerank = 3;
            end
            
            L  = ft_compute_leadfield(eyes, sens, vol, 'reducerank', reducerank);
            
            [sel1, sel2] = spm_match_str(sconf.label, D.chanlabels(chanind));
            
            sconf.coeff(sel1, :) = spm_cond_units(L(sel2, :));
            sconf.bad(sel1, :) = 0;
        end       
    case 'BESA'       
        sconf = spm_eeg_read_bsa(S.conffile );
    case 'SVD'       
        cl = S.conditions;
        svdinput = [];
        for i = 1:numel(cl)
            tmp      = D(D.indchantype('MEEG', 'GOOD'), D.indsample(1e-3*S.timewin(1)):D.indsample(1e-3*S.timewin(2)), D.indtrial(cl{i}));
            svdinput = [svdinput, reshape(tmp, size(tmp, 1),  [])];
        end
        [U, L, V] = spm_svd(svdinput);
        
        
        if isfield(S, 'threshold')
            temp = zeros(size(svdinput));
            for n = size(V, 2):-1:1
                temp = temp+U(:, n)*L(n,n)*V(:,n)';
                if max(max(abs(temp)))>S.threshold;
                    S.ncomp = min(n+1, size(V, 2));
                    break;
                else
                    S.ncomp = 0;
                end
            end       
        end
        
        if S.ncomp>0
            ncomp = min(S.ncomp, size(U, 2));
            [sel1, sel2] = spm_match_str(D.chanlabels(D.indchantype('MEEG')), D.chanlabels(D.indchantype('MEEG', 'GOOD')));
            sconf = [];
            sconf.label = D.chanlabels(D.indchantype('MEEG'));
            sconf.coeff = nan(length(sconf.label), ncomp);
            sconf.coeff(sel1, :) = U(sel2, 1:ncomp);
            sconf.bad = ones(length(sconf.label), 1);
            sconf.bad(sel1, :) = 0;
        end
    case 'SPMEEG'       
        Ds = spm_eeg_load(S.conffile);
        sconf = getfield(Ds, 'sconfounds');        
    case 'CLEAR'
        D = rmfield(D, 'sconfounds');
end

if ~isempty(sconf)
    D = sconfounds(D, sconf, 'append');
end

% Plot scalp topographies
% ---------------------------------------------------------------------
if any(any(D.sconfounds))
    
    Fgraph = spm_figure('GetWin','Graphics');clf
    
    in = [];
    in.f = Fgraph;
    in.noButtons = 1;
    in.cbar = 0;
    in.plotpos = 0;
    
    [junk, modalities] = modality(D, 1, 1);
    
    conf = getfield(D, 'sconfounds');
    
    nm = numel(modalities);
    nc = size(conf.coeff, 2);
    
    for i = 1:nc
        for j = 1:nm
            in.type = modalities{j};
            
            ind = D.indchantype(modalities{j}, 'GOOD');
            
            [sel1, sel2] = spm_match_str(D.chanlabels(ind), conf.label);
            
            Y = conf.coeff(sel2, i);            
            
            in.max = max(abs(Y));
            in.min = -in.max;
            
            in.ParentAxes = subplot(nc, nm, (i - 1)*nm + j);
            spm_eeg_plotScalpData(Y, D.coor2D(ind) , D.chanlabels(ind), in);
            title(sprintf('%s\ncomponent %.0f', modalities{j}, i));           
        end
    end
end

D = D.history(mfilename, S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName', 'Define spatial confounds: done');
