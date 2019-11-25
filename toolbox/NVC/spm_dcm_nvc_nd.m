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
% $Id: spm_dcm_nvc_nd.m 7713 2019-11-25 16:00:34Z spm $

% DCM for MEG simualtion
%--------------------------------------------------------------------------
rng('default')
spm('defaults','EEG');
model          = DCM.model; % Model specification
model{4}       = DCM.N    ; % Vector of included populations
Uf             = DCM.U    ; % fMRI inputs
DCM            = DCM.MEEG ; % M/EEG DCM

%--------------------------------------------------------------------------
Nc               =  size(DCM.C,1);
Ns               =  Nc;
options.spatial  = 'LFP';
options.model    = 'TFM';
if (strcmp(options.model, 'TFM'))
    num_pop      = 4;
end

%--------------------------------------------------------------------------
pE       = DCM.Ep;
[x1]     = spm_dcm_x_neural(pE,options.model);
f        = DCM.M.f ;
xx       = DCM.x{1,1}(1,:);
x        = spm_unvec(xx,x1);
B        = DCM.B{1,1};

% Integrator
%--------------------------------------------------------------------------
M.IS   = 'spm_gen_erp';
M.G    = 'spm_lx_erp';
M.f    = f;
M.x    = x;
M.pE   = pE;
M.m    = length(B);
M.n    = length(spm_vec(M.x));
M.l    = Nc;
M.ns   = DCM.M.ns;
num_condition = size(DCM.xY.y,2);

%--------------------------------------------------------------------------
dt          = DCM.xU.dt;
pst         = (1:M.ns)*dt;
t           =  pst;
M.ons       = DCM.M.ons;
M.dur       = DCM.M.dur;
U.dt        = dt;
U.X         = DCM.xU.X;
P           = pE;
S           = full(Uf.u);
U_fMRI      =  S ;
time_fMRI   = length(S)*Uf.dt ;
in.u        = feval('spm_erp_u',(1:M.ns)*U.dt,P,M);

% Simulation of post-synaptic signals
%--------------------------------------------------------------------------
if (strcmp(model{1}, 'post'))
    i_reg       = 1: Ns            ;
    index_ss    = 8*(i_reg-1) + 1  ;
    index_sp    = index_ss    + 2  ;
    index_inh   = index_sp    + 2  ;
    index_dp    = index_inh   + 2  ;
    
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
    for i= 1: num_condition
        for j = 1: Ns
            R_SS{i,j}  = rms(SS{i,j}) .* U_fMRI(:,i);
            R_SP{i,j}  = rms(SP{i,j}) .* U_fMRI(:,i);
            R_INH{i,j} = rms(INH{i,j}).* U_fMRI(:,i);
            R_DP{i,j}  = rms(DP{i,j}) .* U_fMRI(:,i);
        end
    end
    for  j =1:Ns
        BSS = zeros(size(R_SS{1,1}));
        ASS=BSS ; ASP=BSS ;BSP=BSS ;AINH =BSS ; BINH = BSS ; ADP=BSS ; BDP=BSS ;
        for  i = 1:num_condition
            ASS(:)  =  R_SS{i,j};
            BSS(:)  =  BSS(:)+ASS(:);
            ASP(:)  =  R_SP{i,j};
            BSP(:)  =  BSP(:)+ASP(:);
            AINH(:) =  R_INH{i,j};
            BINH(:) =  BINH(:)+AINH(:);
            ADP(:)  =  R_DP{i,j};
            BDP(:)  =  BDP(:)+ADP(:);
        end
        neuronal_drive.input{j,1} = model{1,4}(1).*BSS;
        neuronal_drive.input{j,2} = model{1,4}(2).*BSP;
        neuronal_drive.input{j,3} = model{1,4}(3).*BINH;
        neuronal_drive.input{j,4} = model{1,4}(4).*BDP;
    end
    neuronal_drive.num = 4;
end

% Simulation of pre-synaptic signals
%--------------------------------------------------------------------------
if (strcmp(model{1}, 'pre'))
    sig = {};
    pq =[]; Q = {};
    for i = 1 : num_condition
        P.xc     = i;
        Q{end+1} = spm_gen_par(P,M,U);
        for j = 1 : size(pst,2)-1
            current_state = DCM.x{i,1}(j,:); %This is results for condition 1
            input = in.u(j);
            [u] = spm_fx_cmc_tfm_gen(current_state,input,Q{1,i},M,model);
            pq(:,:,j) = u;
        end
        sig{i}= pq;
    end
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
    
    RMS_tap_sig = {};
    for i = 1 : num_condition
        for region = 1:Ns
            for pop = 1:num_pop
                RMS_tap_sig{1,i}(region,pop) = rms(sig{1,i}(region,pop,:));
            end
        end
    end
    Rep_sig ={};
    for i = 1: num_condition
        for region = 1: Ns
            for pop = 1:num_pop
                Rep_sig{1,i}(region,pop,: )  = RMS_tap_sig{1,i}(region,pop).*full(U_fMRI(:,i));
            end
        end
    end
    neuronal_drive =[];
    Dr = zeros(size(Rep_sig{1,1}));
    BS = zeros(size(Rep_sig{1,1}));
    for  region =1:Ns
        for pop = 1:num_pop
            for  i = 1:num_condition
                Dr(region,pop,:) =  Rep_sig{1,i}(region,pop,:);
                BS(region,pop,:) =  BS(region,pop,:)+ Dr(region,pop,:);
                Dr=[];
            end
        end
        neuronal_drive.input = BS;
    end
    neuronal_drive.num = 1;
end

% Simulation of decomposed pre-synaptic signals
%--------------------------------------------------------------------------
if (strcmp(model{1}, 'de'))
    sig = {};
    pq =[]; Q = {}; ux = {} ; vx = {}; wx = {};pe = []; pih = [];
    for i = 1 : num_condition
        P.xc     = i;
        Q{end+1} = spm_gen_par(P,M,U);
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
    
    neuronal_drive ={};
    Dr1 = zeros(size(Rep_sig{1,1}));
    Dr2 = Dr1 ;  Dr3 = Dr1 ;  BS1 = Dr1 ; BS2 = Dr1 ;BS3 = Dr1 ;
    for  region =1:Ns
        for pop = 1:num_pop
            for  i = 1:num_condition
                Dr1(region,pop,:)     =  Rep_sig{1,i}(region,pop,:);
                BS1(region,pop,:)     =  BS1(region,pop,:)+ Dr1(region,pop,:);
                Dr1 = [];
                
                Dr2(region,pop,:)     =  Rep_sig{2,i}(region,pop,:);
                BS2(region,pop,:)     =  BS2(region,pop,:)+ Dr2(region,pop,:);
                Dr2 = [];
                
                if (strcmp(model{3}, 'ext'))
                    Dr3(region,pop,:) =  Rep_sig{3,i}(region,pop,:);
                    BS3(region,pop,:) =  BS3(region,pop,:)+ Dr3(region,pop,:);
                    Dr3 = [];
                end
            end
        end
        neuronal_drive.input{1,1}      = BS1;
        neuronal_drive.input{1,2}      = BS1;
        if (strcmp(model{3}, 'ext'))
            neuronal_drive.input{1,3}  = BS3;
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
