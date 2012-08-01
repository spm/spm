function [y,w] = spm_lfp_mtf_sample(P,M,U)
% Spectral response of a NMM (transfer function x noise spectrum)
% FORMAT [G,w] = spm_lfp_mtf_sample(P,M,U)
%
% P - parameters
% M - neural mass model structure
% U - trial-specific effects
%
% G - {G(N,nc,nc}} - cross-spectral density for nc channels {trials}
%                  - for N frequencies in M.Hz [default 1:64Hz]
% w - frequencies
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_lfp_mtf_sample.m 4820 2012-08-01 12:20:00Z guillaume $


% compute log-spectral density
%==========================================================================

% frequencies of interest
%--------------------------------------------------------------------------

% try
%     dt = 1/(2*round(M.Hz(end)));
%     N  = 1/dt;
%     If = round(linspace(M.Hz(1),M.Hz(end),length(M.Hz)));
% catch
%     N  = 128;
%     dt = 1/N;
%     If = 1:N/2;
% end
f    = M.Hz';
w = f*2*pi;


% exogenous (neuronal) inputs
%--------------------------------------------------------------------------
M.u  = sparse(M.m,1);

% solve for fixed point (i.e., 64ms burn in)
%--------------------------------------------------------------------------
S    = M;
S.g  = {};

% initialise states
%-----------------------------------
model = 'NMM';
[x M_int] =  spm_dcm_x_neural_NMDA(P,model);
M.x = x;
   


% channel noise (specific and non-specific)
%--------------------------------------------------------------------------
%Gs   = (exp(P.c(1))*f.^(-1) + exp(P.d(1)))/8;
%Gn   = (exp(P.c(1))*f.^(-1) + exp(P.d(1)))/16;
%Gn   = (exp(P.d(1)))/16;



% cycle over trials
%--------------------------------------------------------------------------
for  c = 1:size(M.X_design,1)
    
    %% spectrum of innovations (Gu)
    %% Independent input for both trials
     Gu   = (exp(P.a(c))*f.^(-1) + exp(P.b(c)));  % spectral density of (AR) input + IID input
  
    
   for i = 1:size(M.X_design,2)
       %% Comment for "first pass" analysis
       
       P.scale_gi= P.scale_gi  +  M.X_design(c,i)*P.B_gi{i};      % intrinsic connections
       P.scale_geP1 = P.scale_geP1 +  M.X_design(c,i)*P.B_geP1{i};
       P.scale_geP2 = P.scale_geP2 +  M.X_design(c,i)*P.B_geP2{i};
       P.scale_geS = P.scale_geS +  M.X_design(c,i)*P.B_geS{i};
       P.scale_geS2 = P.scale_geS2 +  M.X_design(c,i)*P.B_geS2{i};
       P.scale_NMDA = P.scale_NMDA +  M.X_design(c,i)*P.B_NMDA{i};
   
       
  end
    
    M_int = M;
    dt    = 2;
    t     = [1:dt:512]';
    N     = length(t);
    U.dt  = dt/1000;
    U.u   = exp(P.scale_u)*ones(256,1);
    M_int.ns  = 512;
    M_int.f = 'spm_fx_mfm_scale_gaba_NMDA4_dcm';
    M_int.g = [];
    YFM   = spm_int_L(P,M_int,U);
    %%%------------------------------------------------------------ POLE ZERO REP
    % delta
    period{1} = [0:dt/1000:1/2];
    % theta
    period{2} = [0:dt/1000:1/6];
    % alpha
    period{3} = [0:dt/1000:1/12];
    % beta
    period{4} = [0:dt/1000:1/24];
    % gamma
    period{5} = [0:dt/1000:1/46];
    % high gamma
    period{6} = [0:dt/1000:1/80];

    
    for possible = 1:6 %% delta to high gamma
  
        tmp = find(period{possible} > period{possible}(end)/8);
        pyram_samples_spaced(possible,:) = diff(YFM(:,3));
        tmp2 = pyram_samples_spaced(possible,[255:(-tmp(1)):(256 -tmp(1)*7)]);
   
        diff_samples_spaced(possible) = length(find(sign(tmp2)==1));
        samples_spaced(possible,:,:) = YFM([256:(-tmp(1)):(256 -tmp(1)*7)],:);

    end
    

    %% Select evenly spaced samples
    if find(diff_samples_spaced == 4);
        use_set = find(diff_samples_spaced == 4);
    elseif find(diff_samples_spaced == 3);
        use_set = find(diff_samples_spaced == 3);
    else
        use_set = 1;
    end
    x_set = squeeze(samples_spaced(use_set(1),:,:));
   
    
    for sample = 1:8
        
        x_use = squeeze(x_set(sample,:));
        x_jac = spm_unvec(x_use,x);
        
        %%%------------------------------------------------------------ POLE ZERO REP
        M.x = x_jac;
        M.n = 1;
        M.m = 1;
        
        [M0,M1,L] = spm_bireduce(M,P);
        
        
        try
            L = M.U'*L;
        end
        
        
        nc = 1;
        nu = 1;
        A_uc = full(M0(2:end,2:end));
        B_uc = zeros(12,1);
        B_uc(1)=1*32;
        C_uc = zeros(1,12);
        C_uc(3)=1;
        D_uc = 0;

        G     = zeros(8,length(w),nc,nc);
        for i = 1:nc
            for j = 1:nc    
                for k =1:nu
                    [b,a]    =  spm_ss2tf(A_uc,B_uc,C_uc,D_uc);
                    %%%%--------------------------------------------------------   II POWER REP
                    Si2                     =  spm_freqs(b,a,w);
                    Gij                    =  Si2.*conj(Si2);
                    Gij         = (abs(Gij).*Gu)';
                    G(sample,:,i,j) = G(sample,:,i,j) + Gij;
                    G(sample,:,j,i) = G(sample,:,j,i) + Gij;
                    
                end
                %cross-spectral density from channel noise
                %--------------------------------------------------------------
                               G(sample,:,i,j) = G(sample,:,i,j) + Gn';            % common noise
                               save 'G' 'G'
                                save 'Gn' 'Gn'
                                if i == j
                
                                    G(sample,:,i,i) = G(sample,:,i,i) + Gs';       % and channel specific
                                else
                                    % fill in lower half of CSD matrix
                                    %------------------------------------------------------
                                    G(sample,:,j,i) = G(sample,:,i,j);
                end
                
                
            end
        end
    end

    % save frequencies of interest
    % ----------------------------------------------------------------------
    G = squeeze(mean(G));
    y{c} = G;
  end

