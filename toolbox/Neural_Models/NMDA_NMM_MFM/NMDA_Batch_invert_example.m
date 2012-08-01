clear all

       
    %%% Loading theta priors and refitting 
    load('DCM_batch_Memory_trials_1_3_iX_mixedfxfx1') 
    pE = DCM.M.pE; 
    pC = DCM.M.pC; 
    %% Otherwise use: spm_nmm_priors_NMDA.m
    xY = DCM.xY; %% two trial data
    clear DCM
    
    DCM.name =  'TEST_DCM_1';
    DCM.xY   = xY;    
    DCM.A{1} = 0;
    DCM.A{2} = 0;
    DCM.A{3} = 0;
    DCM.B{1} = 1;
    DCM.C    = 1;
    
    DCM.M.Hz = [2:16];                   %% For Model
    DCM.options.Fdcm = [2 16];           %% For Data Extraction
    DCM.options.model = 'NMM';           %% Model to Use
   
    DCM.options.Tdcm = [800 3200];   
    DCM.options.trials = [2 4];
    DCM.M.X_design = [0 1]';
    DCM.M.U = [1];
    DCM.xU.name = {'DA'};
    DCM.options.D = 1;
    DCM.options.han = 0;
    
    DCM.Sname = {'rIFG'};
    DCM.options.spatial = 'LFP'  
    DCM  = spm_dcm_erp_dipfit(DCM, 1);
    Ns   = size(DCM.C,1);                                   % number of sources
    Nc   = DCM.M.dipfit.Nc;
    
    DCM.M.dipfit.model = 'NMM';    
    [x,f]    = spm_dcm_x_neural_NMDA(pE,'NMM');
    
    
    DCM.M.IS = 'spm_lfp_mtf_sample';
    DCM.M.FS = 'spm_lfp_sqrt';
    DCM.M.g  = 'spm_gx_erp';
    DCM.M.f  = 'spm_fx_mfm_NMDA';
    DCM.M.x  = x;
    DCM.M.n  = length(spm_vec(x));
    DCM.M.pE = pE;
    DCM.M.pC = pC;
    DCM.M.pE.B_u = 0;
    DCM.M.m  = 1;
  
    
    [Qp,Cp,Ce,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);
    pE = DCM.M.pE;
    dp  = spm_vec(Qp) - spm_vec(pE);
    Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
    
    % predictions and error (source space)
    %--------------------------------------------------------------------------
    Hc  = feval(DCM.M.IS,Qp,DCM.M,DCM.xU);                   % prediction
    Ec  = spm_unvec(spm_vec(DCM.xY.y) - spm_vec(Hc),Hc);     % prediction error
    
    % store estimates in DCM
    %--------------------------------------------------------------------------
    DCM.Ep = Qp;                   % conditional expectation
    DCM.Cp = Cp;                   % conditional covariance
    DCM.Pp = Pp;                   % conditional probability
    DCM.Hc = Hc;                   % conditional responses (y), channel space
    DCM.Rc = Ec;                   % conditional residuals (y), channel space
    DCM.Ce = Ce;                   % ReML error covariance
    DCM.F  = F;                    % Laplace log evidence
    
 
  save(DCM.name, 'DCM');
     
  
  