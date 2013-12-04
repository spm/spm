function D = spm_eeg_spatial_confounds(S)
% This function defines spatial confounds and adds them to MEEG dataset.
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_spatial_confounds.m 5775 2013-12-04 13:03:55Z vladimir $


SVNrev = '$Rev: 5775 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Define spatial confounds');

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

if ~isfield(S, 'method')
    S.method = spm_input('How to define?','+1', 'm', 'SVD|SPMEEG|Eyes|BESA|Clear', strvcat('SVD', 'SPMEEG', 'Eyes', 'BESA', 'Clear'));
end

switch upper(S.method)
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
        
        eyes = gifti(struct('pnt', [-34 53 -38; 34 53 -38]));
        
        
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
            
            L  = ft_compute_leadfield(eyes.vertices, sens, vol, 'reducerank', reducerank);
            
            [sel1, sel2] = spm_match_str(sconf.label, D.chanlabels(chanind));
            
            sconf.coeff(sel1, :) = spm_cond_units(L(sel2, :));
            sconf.bad(sel1, :) = 0;
        end
        
        D = sconfounds(D, sconf);
    case 'BESA'
        if ~isfield(S, 'conffile')
            S.conffile = spm_select(1, '\.bsa$', 'Select BESA *.bsa file');
        end
        sconf = spm_eeg_read_bsa(S.conffile );
        D = sconfounds(D, sconf);
    case 'SVD'
        if ~isfield(S, 'timewin')
            S.timewin  = spm_input('PST window (ms)', '+1', 'r', '', 2)/1000;
        end
        cl = D.condlist;
        svdinput = [];
        for i = 1:numel(cl)
            svdinput = [svdinput mean(D(D.indchantype('MEEG', 'GOOD'), D.indsample(S.timewin(1)):D.indsample(S.timewin(2)), D.indtrial(cl{i})), 3)];
        end
        [U, L, V] = spm_svd(svdinput);
        
        
        if isfield(S, 'svdthresh')
            temp = zeros(size(svdinput));
            for n = size(V, 2):-1:1
                temp = temp+U(:, n)*L(n,n)*V(:,n)';
                if max(max(temp))>S.svdthresh;
                    S.ncomp = min(n+1, size(V, 2));
                    break;
                else
                    S.ncomp = 0;
                end
            end
        elseif ~isfield(S, 'ncomp')
            S.ncomp = spm_input('How many components?', '+1', 'n', '1', 1);
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
            D = sconfounds(D, sconf);
        end
    case 'SPMEEG'
        if ~isfield(S, 'conffile')
            S.conffile =  spm_select(1, 'mat', 'Select M/EEG mat file with sconfounds');
        end
        Ds = spm_eeg_load(S.conffile);
        sconf = getfield(Ds, 'sconfounds');
        D = sconfounds(D, sconf);
    case 'CLEAR'
        D = rmfield(D, 'sconfounds');
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
