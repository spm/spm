function [D] = spm_eeg_inv_results(D)
% contrast of evoked responses and power for an MEG-EEG model
% FORMAT [D] = spm_eeg_inv_results(D)
% Requires:
%
%     D.inv{i}.contrast.woi   - (n x 2) time (ms) window[s] of interest
%     D.inv{i}.contrast.fboi  - frequency window of interest
%     D.inv{i}.contrast.type  - 'evoked' or 'induced'
%
% this routine will create a contrast for each trial type and will compute
% induced responses in terms of power (over trials) if requested; otherwise
% the power in D.inv{i}.contrast.GW corresponds to the evoked power
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_inv_results.m 3976 2010-07-08 14:12:31Z karl $

% SPM data structure
%==========================================================================
try
    model = D.inv{D.val};
catch
    model = D.inv{end};
    D.val = numel(D.inv);
end

% defaults
%--------------------------------------------------------------------------
try, woi  = model.contrast.woi;  catch, woi  = model.inverse.woi; end
try, foi  = model.contrast.fboi; catch, foi  = [];                end
try, type = model.contrast.type; catch, type = 'evoked';          end


% Ensure contrast woi is within inversion woi
%--------------------------------------------------------------------------
woi(:,1) = max(woi(:,1),model.inverse.woi(1));
woi(:,2) = min(woi(:,2),model.inverse.woi(2));

if ~any(foi)
    foi = [];
end

fprintf('\ncomputing contrast - please wait\n')

% inversion parameters
%--------------------------------------------------------------------------
J    = model.inverse.J;                        % Trial average MAP estimate
T    = model.inverse.T;                        % temporal projector
U    = model.inverse.U;                        % spatial  projector[s]
Is   = model.inverse.Is;                       % Indices of ARD vertices
Ic   = model.inverse.Ic;                       % Indices of channels
It   = model.inverse.It;                       % Indices of time bins
pst  = model.inverse.pst;                      % peristimulus time (ms)
Nd   = model.inverse.Nd;                       % number of mesh dipoles
Nb   = size(T,1);                              % number of time bins
Nc   = size(U,1);                              % number of channels
Nw   = size(woi,1);                            % number of contrast windows
Nj   = numel(J);                               % number of conditions

try
    scale = model.inverse.scale;               % Trial average MAP estimate
catch
    scale = 1;
end


% time-frequency contrast
%==========================================================================
model.contrast.W  = {};
model.contrast.JW = {};
model.contrast.GW = {};

for w = 1:Nw;
    
    % get [Gaussian] time window
    %--------------------------------------------------------------------------
    fwhm = max(diff(woi(w,:)),8);
    t    = exp(-4*log(2)*(pst(:) - mean(woi(w,:))).^2/(fwhm^2));
    t    = t/sum(t);
    
    
    % get frequency space and put PST subspace into contrast (W -> T*T'*W)
    %--------------------------------------------------------------------------
    if ~isempty(foi)
        wt = 2*pi*pst(:)/1000;
        W  = [];
        for f = foi(1):foi(end)
            W = [W sin(f*wt) cos(f*wt)];
        end
        W  = diag(t)*W;
        W  = spm_svd(W,1);
    else
        W  = t(:);
    end
    TW     = T'*W;
    TTW    = T*TW;
    
    % MAP projector and conditional covariance
    %==========================================================================
    qC     = model.inverse.qC*trace(TTW'*model.inverse.qV*TTW);
    qC     = max(qC,0);
    
    
    % cycle over trial types
    %==========================================================================
    try
        trial = model.inverse.trials;
    catch
        trial = D.condlist;
    end
    for i = 1:Nj
        
        % induced or evoked
        %------------------------------------------------------------------
        iw     = (w - 1)*Nj + i;
        CW{iw} = W;
        
        switch(type)
            
            % energy of conditional mean
            %--------------------------------------------------------------
            case{'evoked'}
              
                JW{iw} = J{i}*TW(:,1);
                GW{iw} = sum((J{i}*TW).^2,2) + qC;
                
            % mean energy over trials
            %--------------------------------------------------------------
            case{'induced'}
                
                JW{iw} = sparse(0);
                JWWJ   = sparse(0);
                
                c = D.pickconditions(trial{i});
                
                % conditional expectation of contrast (J*W) and its energy
                %----------------------------------------------------------
                Nt    = length(c);
                spm_progress_bar('Init',Nt,sprintf('condition %d',i),'trials');
                for j = 1:Nt
                    if ~strcmp(D.modality(1,1), 'Multimodal')
                        
                        % unimodal data
                        %--------------------------------------------------
                        Y  = D(Ic{1},It,c(j))*TTW;
                        Y  = U{1}*Y*scale;
                        
                    else
                        
                        % multimodal data
                        %--------------------------------------------------
                        for k = 1:length(U)
                            Y       = D(Ic{k},It,c(j))*TTW;
                            UY{k,1} = U{k}*Y*scale(k);
                        end
                        Y = spm_cat(UY);
                    end
                    
                    MYW    = model.inverse.M*Y;
                    JW{iw} = JW{iw} + MYW(:,1);
                    JWWJ   = JWWJ   + sum(MYW.^2,2);
                    spm_progress_bar('Set',j)
                    
                end
                spm_progress_bar('Clear')
                
                
                % conditional expectation of total energy (source space GW)
                %----------------------------------------------------------
                JW{iw} = JW{iw}/Nt;
                GW{iw} = JWWJ/Nt + qC;
                
            case 'trials'
                
                JW{iw} = {};
                JWWJ   = {};
                
                c = D.pickconditions(trial{i});
                
                % conditional expectation of contrast (J*W) and its energy
                %----------------------------------------------------------
                Nt    = length(c);
                spm_progress_bar('Init',Nt,sprintf('condition %d',i),'trials');
                for j = 1:Nt
                    if ~strcmp(D.modality(1,1), 'Multimodal')
                        
                        % unimodal data
                        %--------------------------------------------------
                        Y     = D(Ic{1},It,c(j))*TTW;
                        Y     = U{1}*Y*scale;
                        
                    else
                        
                        % multimodal data
                        %--------------------------------------------------
                        for k = 1:length(U)
                            Y       = D(Ic{k},It,c(j))*TTW;
                            UY{k,1} = U{k}*Y*scale(k);
                        end
                        Y = spm_cat(UY);
                    end
                    
                    MYW       = model.inverse.M*Y;
                    JW{iw}{j} = MYW(:,1);
                    GW{iw}{j} = sum(MYW.^2,2) + qC;
                    spm_progress_bar('Set',j)
                end
                spm_progress_bar('Clear')
        end
        
    end
    
    
    % Save results
    %======================================================================
    model.contrast.woi   = woi;
    model.contrast.fboi  = foi;
    
    model.contrast.W  = CW;
    model.contrast.JW = JW;
    model.contrast.GW = GW;
    
end % window

D.inv{D.val}         = model;


% Display
%==========================================================================
spm_eeg_inv_results_display(D);
