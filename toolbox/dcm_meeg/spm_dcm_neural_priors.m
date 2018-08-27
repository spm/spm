function [pE,pC] = spm_dcm_neural_priors(A,B,C,model)
% Prepares the priors on the parameters of neural mass models
% FORMAT [pE,pC] = spm_dcm_neural_priors(A,B,C,'model'))
%
% A,B{m},C  - binary constraints on extrinsic connections for m conditions
% 'model'   - 'ERP','SEP','CMC','LFP','NNM' or 'MFM'
%
% pE - prior expectation - f(x,u,P,M)
%
% synaptic parameters (for NMN and MFM)
%--------------------------------------------------------------------------
%    pE.T - synaptic time constants
%    pE.H - synaptic densities
%    pE.S - activation function parameters
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.A  - extrinsic
%    pE.B  - trial-dependent
%    pE.C  - stimulus input
%
%    pE.D  - delays
%
% stimulus and noise parameters
%--------------------------------------------------------------------------
%    pE.R - onset and dispersion
%    pE.U - endogenous activity
%
% pC - prior (co)variances
%
% Because priors are specified under log normal assumptions, most
% parameters are simply scaling coefficients with a prior expectation
% and variance of one.  After log transform this renders pE = 0 and
% pC = 1;
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_neural_priors.m 7409 2018-08-27 11:39:00Z bernadette $

% For generic schemes one can mix and match different types of sources;
% furthermore, they can have different condition-specific modulation of
% intrinsic connectivity parameters and different, source-specific-
% contribution to the lead field (or electrode gain). Source-specific
% models are specified by a structure array model, For the i-th source:
%
% model(i).source  = 'ECD','CMC',...  % source model
% model(i).B       = [i j k ,...]     % free parameters that have B effects
% model(i).J       = [i j k ,...]     % cardinal states contributing to L
% model(i).K       = [i j k ,...]     % other states contributing to L
%__________________________________________________________________________


% assemble priors for more than one sort of model
%==========================================================================
if isstruct(model)
    
    % check source specificaton
    %--------------------------------------------------------------------------
    if numel(model) ~= size(A{1},1)
        error('please enure source number and DCM.option.model are consistent')
    end
    
    % extrinsic parameters (assume CMC form)
    %----------------------------------------------------------------------
    ext     = {'A','B','C','D','M','N','R'};
    int     = {'T','G','M'};
    [pE,pC] = spm_cmc_priors(A,B,C);
    pE      = rmfield(pE,int);
    pC      = rmfield(pC,int);
    
    % intrinsic parameters
    %----------------------------------------------------------------------
    for i = 1:numel(model)
        [iE,iC] = spm_dcm_neural_priors({0,0,0,0},{},0,model(i).source);
        
        % remove extrinsic parameters
        %------------------------------------------------------------------
        for j = 1:numel(ext)
            if isfield(iE,ext(j))
                iE = rmfield(iE,ext{j});
                iC = rmfield(iC,ext{j});
            end
        end
        
        % add condition-specific intrinsic effects (if specified)
        %------------------------------------------------------------------
        if isfield(model(i),'B') && numel(B)
            iC.B(model(i).B) = 1/16;
            iE.B = spm_zeros(iC.B(model(i).B));
            
            if numel(iE.B) > numel(iE.G)
                iE.G = iE.B;
                iC.G = iE.G + iC.G(1);
            end
            
        end
        pE.int{i} = iE;
        pC.int{i} = iC;
        
        if isfield(model(i),'Zf') && iscell(model(i).Zf)
            
            for k = 1:length(pC.B)
                
               if size(model(i).Zf,1)==1
                    pE.Bf{1,k} = zeros(size(model(i).Zf{1,k}));
                    pC.Bf{1,k} = model(i).Zf{1,k} * 1/8;
               elseif size(model(i).Zf,1)==2
                    pE.Bf{1,k} = zeros(size(model(i).Zf{1,k}));
                    pC.Bf{1,k} = model(i).Zf{1,k} * 1/8;
                    pE.Bf{2,k} = zeros(size(model(i).Zf{1,k}));
                    pC.Bf{2,k} = model(i).Zf{2,k} * 1/8;
               else disp('uncorrect specification Zf')
               end
               if size(model(i).Zb,1)==1
                    pE.Bb{1,k} = zeros(size(model(i).Zb{1,k}));
                    pC.Bb{1,k} = model(i).Zb{1,k} * 1/8;
               elseif size(model(i).Zb,1)==2
                    pE.Bb{1,k} = zeros(size(model(i).Zb{1,k}));
                    pC.Bb{1,k} = model(i).Zb{1,k} * 1/8;
                    pE.Bb{2,k} = zeros(size(model(i).Zb{1,k}));
                    pC.Bb{2,k} = model(i).Zb{2,k} * 1/8;
                   else disp('uncorrect specification Zb')
               end
                pC.B{k} = 0*pC.B{k};
            end
        end
    end
    return
end




% get priors on neural model
%--------------------------------------------------------------------------
switch lower(model)
    
    % linear David et al model (linear in states)
    %======================================================================
    case{'erp','sep'}
        
        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_erp_priors(A,B,C);
        
        
        % linear David et al model (Canonical microcircuit)
        %======================================================================
    case{'cmc'}
        
        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_cmc_priors(A,B,C);
        
        % linear David et al model (Canonical microcircuit with plasticity)
        %======================================================================
    case{'tfm'}
        
        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_tfm_priors(A,B,C);
        
        
        % linear David et al model (linear in states)  - with self-inhibition
        %======================================================================
    case{'lfp'}
        
        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_lfp_priors(A,B,C);
        
        
        % Neural mass model (nonlinear in states)
        %======================================================================
    case{'nmm','mfm'}
        
        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_nmm_priors(A,B,C);
        
        % Neural mass model (nonlinear in states)
        %======================================================================
    case{'nmda'}
        
        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_nmda_priors(A,B,C);
        
        % Canonical neural mass model (nonlinear in states)
        %======================================================================
    case{'cmm'}
        
        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_cmm_priors(A,B,C);
        
        % Canonical neural mass model with NMDA channels(nonlinear in states)
        %======================================================================
    case{'cmm_nmda'}
        
        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_cmm_NMDA_priors(A,B,C);
        
        
        % Neural field model (linear in states)
        %======================================================================
    case{'nfm'}
        
        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_nfm_priors(A,B,C);
        
        
        % Neural field model (linear in states)
        %======================================================================
    case{'bgt'}
        
        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_bgt_priors; % no (A,B,C);
        
        
        % Bhatt et al model (motor microcircuit)
        %======================================================================
    case{'mmc'}
        
        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_mmc_priors(A,B,C);
        
        
        % Null model - of Jacabian (linear in states)
        %======================================================================
    case{'null'}
        
        % prior moments on parameters
        %------------------------------------------------------------------
        [pE,pC] = spm_null_priors(A,B,C);
        
        
    otherwise
        warndlg('Unknown model')
end
