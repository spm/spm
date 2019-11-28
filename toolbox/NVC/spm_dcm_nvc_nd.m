function neuronal_drive = spm_dcm_nvc_nd(DCM)
% Generate neuronal drive signals for multimodal DCM for fMRI and M/EEG
% FORMAT neuronal_drive = spm_dcm_nvc_nd(DCM)
%
% Inputs:
% -------------------------------------------------------------------------
% DCM              -  (unestimated multimodal) DCM for fMRI and MEG.
%                     see spm_dcm_nvc_specify.m
%
% Evaluates:
% -------------------------------------------------------------------------
% neuronal_drive   -  neural_drive signals.
%__________________________________________________________________________
% Jafarian, A., Litvak, V., Cagnan, H., Friston, K.J. and Zeidman, P., 2019.
% Neurovascular coupling: insights from multi-modal dynamic causal modelling
% of fMRI and MEG. arXiv preprint arXiv:1903.07478.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Amirhossein Jafarian
% $Id: spm_dcm_nvc_nd.m 7731 2019-11-28 17:53:00Z peter $

% DCM for MEG simualtion
%--------------------------------------------------------------------------
rng('default')
model          = DCM.model; % Model specification
model{4}       = DCM.N    ; % Vector of included populations
Uf             = DCM.U    ; % fMRI inputs
DCM            = DCM.MEEG ; % M/EEG DCM

%--------------------------------------------------------------------------
Ns               =  size(DCM.C,1); % sources
options.spatial  = 'LFP';
options.model    = 'TFM';
if (strcmp(options.model, 'TFM'))
    num_pop      = 4;
end

% Set M/EEG estimated states from condition 1, time 1, as starting values
%--------------------------------------------------------------------------
pE       = DCM.Ep;                              % Parameter structure
[x1]     = spm_dcm_x_neural(pE,options.model);  % Initial states (regions x states)
f        = DCM.M.f ;                            % Neural function
xx       = DCM.x{1,1}(1,:);                     % Predicted activity DCM.x{condition}(sample,state)
x        = spm_unvec(xx,x1);                    % Predicted activity (regions x states)
B        = DCM.B{1,1};                          % Modulatory inputs (forward?)

% Integrator
%--------------------------------------------------------------------------
M.IS   = 'spm_gen_erp';
M.G    = 'spm_lx_erp';
M.f    = f;
M.x    = x;
M.pE   = pE;
M.m    = length(B);
M.n    = length(spm_vec(M.x));
M.l    = Ns;
M.ns   = DCM.M.ns;
num_condition = size(DCM.xY.y,2);

%--------------------------------------------------------------------------
dt          = DCM.xU.dt;                            % EEG sampling rate
pst         = (1:M.ns)*dt;                          % Peri-stimulus time
M.ons       = DCM.M.ons;                            % EEG onset
M.dur       = DCM.M.dur;                            % EEG duration
U.dt        = dt;
U.X         = DCM.xU.X;
P           = pE;
U_fMRI      = full(Uf.u);                           % fMRI inputs
in.u        = feval('spm_erp_u',(1:M.ns)*U.dt,P,M); % EEG input function

% Simulation of post-synaptic signals
%--------------------------------------------------------------------------
if (strcmp(model{1}, 'post'))
    % Compute indices within the state vector (see spm_fx_cmc_tfm_gen)
    i_reg       = 1: Ns            ;
    index_ss    = 8*(i_reg-1) + 1  ;
    index_sp    = index_ss    + 2  ;
    index_inh   = index_sp    + 2  ;
    index_dp    = index_inh   + 2  ;    
    % Apply Hanning window
    n = size(pst,2);
    tap  = hanning(n, 'symmetric');    
    for i= 1: num_condition
        for j = 1: Ns
            for time = 1:size(pst,2)
                SS{i,j}(time)    =(DCM.x{i,1}(time,index_ss(j)))' *tap(time);
                SP{i,j}(time)    =(DCM.x{i,1}(time,index_sp(j)))' *tap(time);
                INH{i,j}(time)   =(DCM.x{i,1}(time,index_inh(j)))'*tap(time);
                DP{i,j}(time)    =(DCM.x{i,1}(time,index_dp(j)))' *tap(time);
            end
        end
    end
    % Compute root mean squared and scale fMRI inputs
    for i= 1:num_condition
        for j = 1: Ns
            R_SS{i,j}  = rms(SS{i,j}) .* U_fMRI(:,i);
            R_SP{i,j}  = rms(SP{i,j}) .* U_fMRI(:,i);
            R_INH{i,j} = rms(INH{i,j}).* U_fMRI(:,i);
            R_DP{i,j}  = rms(DP{i,j}) .* U_fMRI(:,i);
        end
    end
    % Sum over conditions
    pincluded  = model{4};
    for j = 1:Ns
        BSS = zeros(size(R_SS{1,1}));
        BSP=BSS ; BINH = BSS ; BDP=BSS ;
        for  i = 1:num_condition
            BSS(:)  =  BSS(:)+R_SS{i,j};
            BSP(:)  =  BSP(:)+R_SP{i,j};
            BINH(:) =  BINH(:)+R_INH{i,j};
            BDP(:)  =  BDP(:)+R_DP{i,j};
        end
        neuronal_drive.input{j,1} = pincluded(1).*BSS;
        neuronal_drive.input{j,2} = pincluded(2).*BSP;
        neuronal_drive.input{j,3} = pincluded(3).*BINH;
        neuronal_drive.input{j,4} = pincluded(4).*BDP;
    end
    neuronal_drive.num = 4;
end

% Simulation of pre-synaptic signals
%--------------------------------------------------------------------------
if (strcmp(model{1}, 'pre'))
    % Build matrix of signals [regions x populations x time] {conditions}
    sig = {};
    pq =[]; Q = {};
    for i = 1 : num_condition
        P.xc     = i;
        Q{end+1} = spm_gen_par(P,U);
        for j = 1 : size(pst,2)-1
            % Get hidden state (see spm_fx_cmc_tfm_gen)
            current_state = DCM.x{i,1}(j,:);
            input = in.u(j);
            % Get presynaptic responses
            [u] = spm_fx_cmc_tfm_gen(current_state,input,Q{1,i},M,model);
            pq(:,:,j) = u;
        end
        sig{i}= pq;
    end
    % Apply Hanning window
    n = size(pst,2);
    tap  = hanning(n, 'symmetric');
    tap_sig={};
    for i = 1 : num_condition
        for region = 1:Ns
            for pop = 1:num_pop
                for time = 1: size(pst,2)-1
                    tap_sig{1,i}(region,pop,time) = sig{1,i}(region,pop,time)* tap(time);
                end
            end
        end
    end
    % Convert to root mean squared (collapsing over time)
    RMS_tap_sig = {};
    for i = 1 : num_condition
        for region = 1:Ns
            for pop = 1:num_pop
                RMS_tap_sig{1,i}(region,pop) = rms(sig{1,i}(region,pop,:));
            end
        end
    end
    % Use RMS from EEG to scale fMRI input
    Rep_sig ={};
    for i = 1: num_condition
        for region = 1: Ns
            for pop = 1:num_pop
                Rep_sig{1,i}(region,pop,: )  = RMS_tap_sig{1,i}(region,pop).*full(U_fMRI(:,i));
            end
        end
    end
    % Sum over conditions to create neuronal drive    
    Dr = zeros(Ns,num_pop);
    for c = 1:num_condition
        Dr = Dr + Rep_sig{c};
    end
    neuronal_drive = [];
    neuronal_drive.input = Dr;
    neuronal_drive.num   = 1;
end

% Simulation of decomposed pre-synaptic signals
%--------------------------------------------------------------------------
if (strcmp(model{1}, 'de'))
    % Build matrix [regions x populations x time] {input type x conditions}
    sig = {};
    Q = {}; pe = []; pih = []; pex = [];
    for i = 1 : num_condition
        P.xc     = i;
        Q{end+1} = spm_gen_par(P,U);
        for j = 1 : size(pst,2)-1
            current_state  = DCM.x{i,1}(j,:);
            input          = in.u(j);
            [ux, vx, wx]   = spm_fx_cmc_tfm_gen(current_state,input,Q{1,i},M,model);
            pe(:,:,j)      =  ux;
            pih(:,:,j)     =  vx;
            if (strcmp(model{3}, 'ext'))
                pex(:,:,j) = wx;
            end
        end
        sig{1,i} = pe;
        sig{2,i} = pih;
        if (strcmp(model{3}, 'ext'))
            sig{3,i}= pex;
        end
    end
    % Apply Hanning window
    n       = size(pst,2);
    tap     = hanning(n, 'symmetric');
    tap_sig ={};
    for i = 1 : num_condition
        for region = 1:Ns
            for pop = 1:num_pop
                for time = 1: size(pst,2)-1
                    tap_sig{1,i}(region,pop, time)     = sig{1,i}(region,pop,time)* tap(time);
                    tap_sig{2,i}(region,pop, time)     = sig{2,i}(region,pop,time)* tap(time);
                    if (strcmp(model{3}, 'ext'))
                        tap_sig{3,i}(region,pop, time) = sig{3,i}(region,pop,time)* tap(time);
                    end
                end
            end
        end
    end
    % Compute RMS
    for i = 1 : num_condition
        for region = 1:Ns
            for pop = 1:num_pop
                Rms_sig{1,i}(region,pop)     = rms(sig{1,i}(region,pop,:));
                Rms_sig{2,i}(region,pop)     = rms(sig{2,i}(region,pop,:));
                if (strcmp(model{3}, 'ext'))
                    Rms_sig{3,i}(region,pop) = rms(sig{3,i}(region,pop,:));
                end
            end
        end
    end
    % Use RMS from EEG to scale fMRI input
    Rep_sig ={};
    for i = 1: num_condition
        for region = 1: Ns
            for pop = 1:num_pop
                Rep_sig{1,i}(region,pop,:)      = Rms_sig{1,i}(region,pop).*full(U_fMRI(:,i));
                Rep_sig{2,i}(region,pop,:)      = Rms_sig{2,i}(region,pop).*full(U_fMRI(:,i));
                if (strcmp(model{3}, 'ext'))
                    Rep_sig{3,i}(region,pop,:)  = Rms_sig{3,i}(region,pop).*full(U_fMRI(:,i)) ;
                end
            end
        end
    end
    % Sum over conditions
    neuronal_drive ={};
    Dr1 = zeros(size(Rep_sig{1,1}));
    BS1 = Dr1; BS2 = Dr1; BS3 = Dr1 ;
    for region =1:Ns
        for pop = 1:num_pop
            for i = 1:num_condition
                BS1(region,pop,:) = BS1(region,pop,:) + Rep_sig{1,i}(region,pop,:);
                BS2(region,pop,:) = BS2(region,pop,:) + Rep_sig{2,i}(region,pop,:);
               
                if (strcmp(model{3}, 'ext'))
                    BS3(region,pop,:) = BS3(region,pop,:)+ Rep_sig{3,i}(region,pop,:);
                end
            end
        end
        neuronal_drive.input{1,1} = BS1;
        neuronal_drive.input{1,2} = BS2;
        if (strcmp(model{3}, 'ext'))
            neuronal_drive.input{1,3} = BS3;
        end
    end
    if (strcmp(model{3}, 'ext'))
        neuronal_drive.num = 3;
    else
        neuronal_drive.num = 2;
    end
end

end

function y = rms(x)
% Root mean squared value.
    y = sqrt(mean(x .* conj(x)));
end
