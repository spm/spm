function SCKS = spm_SCK(SCKS)
% FORMAT SCKS = spm_SCK(SCKS)
%__________________________________________________________________________
% Square-root Cubature Kalman Filters [2] & Square-root Rauch-Tang-Striebel
% Smoother (SCKF-SCKS [1]).
%==========================================================================
% This function performs joint estimation of the states, input and parameters
% of a model that is described as a stochastic continuous-discrete 
% state-space in terms of nonlinear blind deconvolution. The state equations
% must have the form of ordinary differential equations, where the 
% discretization is performed through local-linearization scheme [3]. 
% Additionally, the parameter noise covariance is estimated online via 
% stochastic Robbins-Monro approximation method [4], and the measurement noise 
% covariance is estimated using a combined variational Bayesian (VB) 
% approach with a nonlinear filter/smoother [5].
%__________________________________________________________________________
%
% SCKS.M  - model structure (based on DEM [6] in SPM8 toolbox)
% SCKS.Y  - response variable, output or data
%__________________________________________________________________________
%
% generative model:
%--------------------------------------------------------------------------
%   M(1).f  = dx/dt = f(x,v,P)    {inline function, string or m-file}
%   M(1).g  = y(t)  = g(x,v,P)    {inline function, string or m-file}
%   
%   M(1).xP = state error covariance matrix
%   M(1).uP = input error variance
%   M(1).wP = parameter error covariance matrix
%
%   M(1).pE = prior expectation of p model-parameters
%   M(1).pC = prior covariances of p model-parameters
%   M(1).ip = parameter indices
%   M(1).cb = constrain on parameters [lower, upper];
%
%   M(1).Q  = precision components on observation noise
%   M(1).V  = fixed precision (input noise)
%   M(1).W  = precision on state noise (approximated by annealing)
%
%   M(i).m  = number of inputs v(i + 1);
%   M(1).n  = number of states x(i);
%   M(1).l  = number of output v(i);
%
%   M(1).Qf      = form of measurement noise cov estimate:
%                  'auto'[default],'min','mean'
%   M(1).E.nN    = number of SCKF-SCKS algorithm iterations
%   M(1).E.Itol  = tolerance value for SCKF-SCKS convergence 
%   M(1).E.nD    = number of integration step between observations
%   M(1).VB.N    = number of VB algorithm iterations
%   M(1).VB.Itol = tolerance value for VB convergence 
%   M(1).VB.l    = VB scaling factor;
%
%   conditional moments of model-states - q(u)
%--------------------------------------------------------------------------
%   qU.x  = Conditional expectation of hidden states (backward estimate)
%   qU.v  = Conditional expectation of input (backward estimate)
%   qU.z  = Conditional prediction error 
%   qU.S  = Conditional covariance: cov(x) (states - backward estimate)
%   qU.C  = Conditional covariance: cov(u) (input - backward estimate)
%
% conditional moments of model-parameters - q(p)
%--------------------------------------------------------------------------
%   qP.P    = Conditional expectation
%   qP.C    = Conditional covariance
%
%      F    = negative log-likelihood
%__________________________________________________________________________
% Copyright (c) Brno University of Technology (2010)...
% Martin Havlicek 05-12-2010
% 
% References:
% [1] Havlicek M et al (2011)
% [2] Arasaratnam, I., Haykin, S. (2009) Cubature Kalman Filters. IEEE
%     Transactions on Automatic Control 54, 1254-1269.
% [3] Jimenez, J.C. (2002) A simple algebraic expression to evaluate the
%     local linearization schemes for stochastic differential equations* 
%     1. Applied Mathematics Letters 15, 775-780.
% [4] Van der Merwe, R., 2004. Sigma-point Kalman filters for probabilistic
%     inference in dynamic state-space models. Ph.D.thesis, Oregon Graduate 
%     Institute of Science and Technology.
% [5] Sarkka, S., Hartikainen, J. (2011?) Extension of VB-AKF to Estimation
%     of Full Covariance and Non-Linear Systems. In Press.
% [6] Friston, K.J., et al. (2008) DEM: a variational treatment of dynamic
%     systems. Neuroimage 41, 849-885.
%__________________________________________________________________________
% Copyright (C) - Martin Havlicek
 
% Martin Havlicek
% $Id: spm_SCK.m 4628 2012-01-27 20:51:41Z karl $
% check model specification
%--------------------------------------------------------------------------
M  = SCKS.M;
M  = spm_DEM_M_set(M);
 
% get integration step dt:
%--------------------------------------------------------------------------
nD = M(1).E.nD;    
dt = 1/nD;            

% INITIALISATION:
% =========================================================================

% interpolate observation according to integration step
%--------------------------------------------------------------------------
y     = SCKS.Y;            % observations
if size(y,1) > size(y,2)   % check the dimensions
    y = y';
end

% interpolate if dt < 1:
%--------------------------------------------------------------------------
y    = spm_interp(y',1/dt)';
if size(y,1) > size(y,2)   % check dimensions again
    y = y';
end
T    = size(y,2);          % number of time points 
 
% initial condition:
%--------------------------------------------------------------------------
x     = M(1).x;            % states
u     = M(2).v;            % input
pE    = spm_vec(M(1).pE);  % all model parameter
ip    = M(1).ip;           % parameter indices to be estimated
theta = pE(ip);            % selected parameters
 
try cb  = M(1).cb;               catch, cb = []; end; % parameter constraints  
try tE  = spm_vec(SCKS.pP.P{1}); catch, tE = []; end; % true parameters for display (if available)
                            
% covariances (square-roots)
%--------------------------------------------------------------------------
sR      = cell(1,T);
[sR{:}] = deal(sparse(real(chol(inv(M(1).V)))*dt));   % observation noise variance
sQ      = sparse(real(chol(inv(M(1).W)))*dt);         % hidden state noise variance
if ~isempty(M(2).v)
    sV  = sparse(real(chol(inv(M(2).V)))*dt);         % input noise variance
else
    sV  = [];
end
                                                           
% process error covariances (square-roots)
%--------------------------------------------------------------------------
Sx = sparse(real(chol(M(1).xP))*dt);
if ~isempty(M(2).v)
    Su = sparse(real(chol(M(1).uP))*dt);
else
    Su = [];
end
if ~isempty(ip)
    Sw = sparse(real(chol(M(1).wP(ip,ip)))*dt);
    sW = sparse(real(chol(M(1).pC(ip,ip)))*dt);  % parameter noise variance
    dv = diag(sW); 
else
    Sw = [];
    sW = [];
end
 
% number of states, inputs and parameters:
%--------------------------------------------------------------------------
nx      = size(Sx,1);         % number of states
nu      = size(sV,1);         % number of states
nw      = size(sW,1);         % number of parameters
no      = size(sR{1},1);      % number of observations
noises  = nx + nu + nw + no;  % number of noise components
 
% concatenate state vector and square-root error covariance:
%--------------------------------------------------------------------------
xc      = [x(:); u(:); theta(:)];
xx      = zeros(nx+nu+nw,T);
xx(:,1) = xc;
Sc      = cell(1,T);
[Sc{:}] = deal(sparse(nx+nu+nw,nx+nu+nw));
Sc{1}   = blkdiag(Sx,Su,Sw);
 
% get vector indices for components of concatenated state vector
%--------------------------------------------------------------------------
xmask   = [ones(1,nx),ones(1,nu)*2,ones(1,nw)*3,ones(1,no)*4];
xind    = find(xmask==1);
uind    = find(xmask==2);
wind    = find(xmask==3);
clear xmask;
 
% setting for VB: observation noise estimation:
%--------------------------------------------------------------------------
if ~isempty(M(1).Q)
    try, iter0  = M(1).VB.N;   catch, iter0 = 3;          end
    try, lambda = M(1).VB.l;   catch, lambda = 1-exp(-2); end
    NU      = 6;
    V       = diag(repmat(1e-4,1,no));
    [sR{:}] = deal(sqrt(1./(NU-no-1).*V));
    B       = sqrt(lambda)*eye(no);
    k       = size(sR{1},1);
    iter    = iter0;
    MSE0    = zeros(no,1);
    RR0     = zeros(no,T);
    VBrun   = [];
else
    iter0   = 1;
    iter    = iter0;
    RR      = [];
    VBrun   = [];
end
 
% Pre-calculate cubature points: 
%--------------------------------------------------------------------------
n          = nx + nu + nw + noises;       % total state vector dimension
nPts       = 2*n;                         % number of cubature points
CubPtArray = sqrt(n)*[eye(n) -eye(n)];    % cubature points array
 
% augment paramter matrix by number of cubature points:
%--------------------------------------------------------------------------
pE   = pE(:,ones(1,nPts));
 
% prepare matrix template for integration by Local linearization scheme:
%--------------------------------------------------------------------------
EXPm = repmat({[ones(nx),2*ones(nx,1);zeros(1,nx+1)]},1,nPts);
EXPm = sparse(blkdiag(EXPm{:}));
xt   = repmat([zeros(1,nx) 1],1,nPts)';
 
 
% =========================================================================
% Iteration scheme:
% =========================================================================
% get maximum number of iterations and tolerance:
%--------------------------------------------------------------------------
try, ItolVB = M(1).VB.Itol;  catch,  ItolVB = 1e-4;      end
try, Itol   = M(1).E.Itol;   catch,  Itol   = 1e-3;      end
try, RUN    = M(1).E.nN;     catch,  RUN    = 32;        end
try, ap     = M(1).E.RM;     catch,  ap     = [1e3 1e6]; end   % Robins-Monro approximation parameters
 
MLdiff0  = 1e-1;
mloglik0 = 0;
ML       = [];
VBrun    = RUN;
t0       = tic;
% =========================================================================
% Iteration loop (until convergence)
% =========================================================================
 
for run = 1:RUN
    t1  = tic;
    mloglik = -log(2.*pi).*(T/dt);
    % =====================================================================
    %   Forward pass:
    % =====================================================================
     for t = 2:T
         
        sQ = diag(diag((1/sqrt(0.9995)-1)*Sc{t-1}(xind,xind)));
        Sa = blkdiag(Sc{t-1},sQ,sV,sW,sR{t-1});
        xa = [xc;zeros(noises,1)];
        Xi =  xa(:,ones(1,nPts)) + Sa*CubPtArray;
        
        %------------------------------------------------------------------
        % PREDICTION STEP:
        %------------------------------------------------------------------
        xPred(uind,:) = Xi(uind,:) + Xi(uind+nx+nu+nw,:);  % add input noise
        xPred(wind,:) = Xi(wind,:) + Xi(wind+nx+nu+nw,:);  % add parameter noise
        
        % parameter constraint:
        %------------------------------------------------------------------
        if ~isempty(cb) && ~isempty(ip)
            xPred(wind,:) = min(cb(:,2*ones(1,nPts)),xPred(wind,:)); % upper constrain
            xPred(wind,:) = max(cb(:,1*ones(1,nPts)),xPred(wind,:)); % lower constrain
        end
        pE(ip,:) = xPred(wind,:);
        
        % propagation of cubature points through nonlinear function:
        %------------------------------------------------------------------
        f             = M(1).f(Xi(xind,:),xPred(uind,:),pE);
        
        % integration by local-linearization scheme:
        %------------------------------------------------------------------
        dfdx          = spm_diff_all(M(1).f,Xi(xind,:),xPred(uind,:),pE,1);
        dx            = expmall(dfdx,f,dt,EXPm)*xt;
        xPred(xind,:) = Xi(xind,:) + reshape(dx(~xt),nx,nPts) + Xi(xind+nx+nu+nw,:);
        
        % mean prediction:
        %------------------------------------------------------------------
        x1      = sum(xPred,2)/nPts;
        X       = (xPred-x1(:,ones(1,nPts)))/sqrt(nPts);
        [foo,S] = qr([X]',0);
        S       = S';
        
        % in the case of VB observation noise estimation:
        %------------------------------------------------------------------
        if ~isempty(M(1).Q) && iter~=1
            NU   = lambda.*(NU-k-1)+k+1;
            V    = B*V*B';
            NU   = NU + 1;
            V0   = V;
        end
        x10    = x1;
        S0     = S;
        yPred0 = M(1).g(xPred(xind,:),xPred(uind,:),pE);
        
        %------------------------------------------------------------------
        % UPDATE STEP:
        %------------------------------------------------------------------
        
        % VB estimation of sR (iteratively)
        %------------------------------------------------------------------
        for it = 1:iter
            
           % VB part - update of square-root measurement noise cov
           %---------------------------------------------------------------
            if ~isempty(M(1).Q) && iter ~= 1
                sR{t} = real(chol(1/(NU - no - 1)*V)); 
            end
            
            % propagate cubature points through observation function:
            %--------------------------------------------------------------
            yPred    = yPred0 + sR{t}*CubPtArray(end-no+1:end,:);
            y1       = sum(yPred,2)/nPts;
            
            Y        = (yPred-y1(:,ones(1,nPts)))/sqrt(nPts);
            [foo,Sy] = qr([Y]',0);
            Sy       = Sy';
            
            Pxy      = X*Y';                % cross covariance
            K        = (Pxy/Sy')/Sy;        % Kalman gain
            resid    =  y(:,t) - y1;        % innovations
            
            % state (input,parameter) estimates:
            %--------------------------------------------------------------
            xc = x10 + K*(resid);
            
            % check parameter constraints: 
            %--------------------------------------------------------------
            if ~isempty(cb) && ~isempty(ip)
                xc(wind) = min(cb(:,2),xc(wind)); % upper constrain
                xc(wind) = max(cb(:,1),xc(wind)); % lower constrain
            end
            % estimate of process error covarinace:
            %--------------------------------------------------------------
            [foo,S] = qr([(X - K*Y)]',0);
            S       = S';
            
            % VB part
            %--------------------------------------------------------------
            if ~isempty(M(1).Q) && iter~=1
                Sa     = blkdiag(S,sQ,sV,sW,sR{t});
                Xi     = repmat([xc;zeros(noises,1)],1,nPts) + Sa*CubPtArray;
                yPreds = M(1).g(Xi(xind,:),Xi(uind,:),pE);  % no additive noise!
                D      = repmat(y(:,t),1,nPts)-yPreds;
                D      = D*D'/nPts;
                V      = V0 + D;
            end
         
        end
        Sc{t}   = S;
        xx(:,t) = xc;
 
        % Maximum log-Likelihood
        %------------------------------------------------------------------
         mloglik = mloglik - log(det(Sy*Sy')) - resid'/(Sy*Sy')*resid;
       
        % Robins-Monro stochastic approximation for of parameters noise cov
        %------------------------------------------------------------------
         if ~isempty(ip)
             subKG = K(wind,:);
             dv    = sqrt((1-1/ap(1))*(dv.^2) + 1/ap(1)*diag(subKG*(subKG*resid*resid')'));
             sW    = diag(dv);
             ap(1) = min(ap(1)+1,ap(2));
         end
         if ~isempty(M(1).Q) && iter~=1
             RR(:,t) = diag(sR{t});
         end
    end
    xxf = xx;
    Sf  = Sc;
    %----------------------------------------------------------------------
    % END of forward pass
    % ---------------------------------------------------------------------
  
    % =====================================================================
    %   Backward pass:
    % =====================================================================
    for t = T-1:-1:1
        % VB part:
        %------------------------------------------------------------------
        if ~isempty(M(1).Q) && iter~=1
            NU    = lambda.*(NU-k-1)+k+1;
            V     = B*V*B';
            NU    = NU + 1;
            V0    = V;
            sR{t} = real(chol(1/(NU-no-1)*V)); 
        end
        
        % Square-root Cubature Rauch-Tung-Striebel smoother
        %------------------------------------------------------------------
        
        % evaluate cubature points:
        %------------------------------------------------------------------
        Sa = blkdiag(Sc{t},sQ,sV,sW,sR{t});
        xa = [xx(:,t);zeros(noises,1)];
        Xi = xa(:,ones(1,nPts)) + Sa*CubPtArray;
      
        xPred(uind,:) = Xi(uind,:) + Xi(uind+nx+nu+nw,:);
        xPred(wind,:) = Xi(wind,:) + Xi(wind+nx+nu+nw,:);
        
        % check parameter constraints:
        %------------------------------------------------------------------
        if ~isempty(cb) && ~isempty(ip)
            xPred(wind,:) = min(cb(:,2*ones(1,nPts)),xPred(wind,:)); % upper constrain
            xPred(wind,:) = max(cb(:,1*ones(1,nPts)),xPred(wind,:)); % lower constrain
        end
        
        pE(ip,:)      = xPred(wind,:);
        % propagate cubature points through nonlinear function:
        %------------------------------------------------------------------
        f             = M(1).f(Xi(xind,:),xPred(uind,:),pE);
        dfdx          = spm_diff_all(M(1).f,Xi(xind,:),xPred(uind,:),pE,1);
        dx            = expmall(dfdx,f,dt,EXPm)*xt;
        xPred(xind,:) = Xi(xind,:) + reshape(dx(~xt),nx,nPts) + Xi(xind+nx+nu+nw,:) ;
    
        x1       = sum(xPred,2)/nPts;
        X        = (xPred-x1(:,ones(1,nPts)))/sqrt(nPts) + eps;
        x01      = xx(:,t);
        X01      = (Xi([xind,uind,wind],:) - repmat(x01,1,nPts))/sqrt(nPts);
        [foo,S]  = qr([X]',0);
        S        = S';
        
        Pxy      = X01*X';      % cross covariance
        K        = (Pxy/S')/S;  % Kalman gain
        
        % smoothed estimate of the states (input, parameters)
        % and process error covariance:
        %------------------------------------------------------------------
        xx(:,t)  = xx(:,t) + K*(xx(:,t+1) - x1);
        [foo,S]  = qr([X01 - K*X, K*Sc{t+1}]',0);
        S        = S';
        Sc{t}    = S;
        
        % check parameter constraint:
        %------------------------------------------------------------------
        if ~isempty(cb) && ~isempty(ip)
            xx(wind,t) = min(cb(:,2),xx(wind,t)); % upper constrain
            xx(wind,t) = max(cb(:,1),xx(wind,t));
        end
        
        % VB part (smoothing):
        %------------------------------------------------------------------
        if ~isempty(M(1).Q) && iter~=1
            Sa     = blkdiag(S,sQ,sV,sW,sR{t});
            Xi     = repmat([xx(:,t);zeros(noises,1)],1,nPts) + Sa*CubPtArray;
            yPreds = M(1).g(Xi(xind,:),Xi(uind,:),pE);  % no additive noise!
            D      = repmat(y(:,t),1,nPts)-yPreds;
            D      = D*D'/nPts;
            V      = V0 + D;
        end
        
        if ~isempty(M(1).Q) && iter~=1
            RR(:,t) = diag(sR{t});
        end
     
    end
    xxb = xx;
    Sb  = Sc;
    %----------------------------------------------------------------------
    % END of backward pass
    %----------------------------------------------------------------------
 
    str{1} = sprintf('SCKS: %i (1:%i)',run,iter);
 
    % iteration condition for measurement noise estimate:
    % iterate until stabilization of sR estimate
    %----------------------------------------------------------------------
    if ~isempty(M(1).Q) && iter0 ~= 1
        
        MSE     = mean((RR -(RR0)).^2,2);
        RR0     = RR;
        MSEdiff = abs(MSE - MSE0);
        MSE0    = MSE;
 
        if MSEdiff < ItolVB  % (till it gets stable)
            switch(lower(M(1).Qf))
                case('all')
                    % take all
                case('mean')
                    sR      = cell(1,T);
                    [sR{:}] = deal(diag(mean(RR,2)));
                case('min')
                    RRs     = sort(RR,2,'descend');
                    sR      = cell(1,T);
                    [sR{:}] = deal(diag(mean(RRs(:,round(T*0.90):end),2)));
                case('auto')
                    dlim    = min(RR,[],2);
                    ulim    = max(RR,[],2);
                    if all(ulim./dlim<4)
                        % take all
                    else
                        RRs     = sort(RR,2,'descend');
                        sR      = cell(1,T);
                        [sR{:}] = deal(diag(mean(RRs(:,round(T*0.90):end),2)));
                    end
            end
            iter0    = 1;
            iter     = iter0;
            mloglik0 = 0;
            VBrun    = run;
        end
    end
    
    
    % log-likelihood difference:
    %----------------------------------------------------------------------
    MLdiff(run) = abs(mloglik-mloglik0);
    ML(run)     = mloglik;
    if mloglik > 0
        mloglik     = mloglik - 5000;
        MLdiff(run) = abs(mloglik-mloglik0);
        ML(run)     = mloglik;
    end
    
    timed  = toc(t1);
    str{2} = sprintf('F:%.4e',ML(end));
    str{3} = sprintf('dF:%.4e',MLdiff(end));
    str{4} = sprintf('(%.2e sec)',timed);
    fprintf('%-16s%-16s%-16s%-16s\n',str{:})
 
    % plot estimates:
    %----------------------------------------------------------------------
    doplotting(M,xxf,xxb,Sf,Sb,ML,T,wind,ip,run,RR,VBrun);
 
    
    % stopping condition:
    %----------------------------------------------------------------------
    if RUN > 1 && (~isempty(ip) || ~isempty(M(1).Q))
        if run == 2
            MLdiff0 = MLdiff(run);
        elseif run > 2
            if MLdiff0 < MLdiff(run),
                MLdiff0 = MLdiff(run);
            end
        end
        if (((MLdiff(run)/MLdiff0)<Itol || run==RUN) || (isempty(ip)&&MLdiff(run)<Itol)) && iter0 == 1,
 
            yy       = M(1).g(xx(xind,:),xx(uind,:),pE(:,1));
            res      = y - yy;
 
            try SCKS = rmfield(SCKS,'qU'); end
            try SCKS = rmfield(SCKS,'qP'); end
            try SCKS = rmfield(SCKS,'qH'); end
            
            % save results in structure:
            %--------------------------------------------------------------
            SCKS.qU.v{1} = yy(:,1:nD:end);
            SCKS.qU.x{1} = xxb(xind,1:nD:end);
            SCKS.qU.v{2} = xxb(uind,1:nD:end);
            SCKS.qU.z{1} = res(:,1:nD:end);
            if ~isempty(ip)
                qP           = SCKS.M(1).pE;
                qP(ip)       = mean(xx(wind,:),2);
                SCKS.qP.P{1} = qP;
                SCKS.qP.P{2} = [];
            end
            SCKS.F       = ML;
            
            for i = 1:nD:T
                j = 1 + (i - 1)/nD;
                SCKS.qU.S{j} = Sb{i}(xind,xind);
                SCKS.qU.C{j} = Sb{i}(uind,uind);
                if ~isempty(ip)
                    SCKS.qP.p{j} = xxb(wind,:);
                    SCKS.qP.c{j} = Sb{i}(wind,wind);
                    SCKS.qP.ip   = ip;
                end
            end
            return
        end
 
        mloglik0 = mloglik;
        xc       = [xx([xind,uind],1); mean(xx(wind,:),2)];
        xx(:,1)  = xc;
 
    else
        
        pE(ip,1) = mean(xx(wind,:),2);
        yy       = M(1).g(xx(xind,:),xx(uind,:),pE(:,1));
        res      = y - yy;
 
        try SCKS = rmfield(SCKS,'qU'); end
        try SCKS = rmfield(SCKS,'qP'); end
        try SCKS = rmfield(SCKS,'qH'); end
        
        % save results in structure:
        %------------------------------------------------------------------
        SCKS.qU.v{1} = yy(:,1:nD:end);
        SCKS.qU.x{1} = xxb(xind,1:nD:end);
        SCKS.qU.v{2} = xxb(uind,1:nD:end);
        SCKS.qU.z{1} = res(:,1:nD:end);
        if ~isempty(ip)
            qP           = SCKS.M(1).pE;
            qP(ip)       = mean(xx(wind,:),2);
            SCKS.qP.P{1} = qP;
            SCKS.qP.P{2} = [];
        end
        SCKS.F       = ML;
        
        for i = 1:nD:T
            j = 1 + (i - 1)/nD;
            SCKS.qU.S{j} = Sb{i}(xind,xind);
            SCKS.qU.C{j} = Sb{i}(uind,uind);
            if ~isempty(ip)
                SCKS.qP.p{j} = xxb(wind,:);
                SCKS.qP.c{j} = Sb{i}(wind,wind);
                SCKS.qP.ip   = ip;
            end
        end
        return
    end
end
 
%==========================================================================
%==========================================================================
 
 
%--------------------------------------------------------------------------
% Plot estimates at each iteration:
%--------------------------------------------------------------------------
function doplotting(M,xxf,xxb,Sf,Sb,ML,T,wind,ip,run,RR,VBrun)
 
% Initialize display:
%--------------------------------------------------------------------------
spm_figure('GetWin','SCKF-SCKS estimates');
set(gcf,'Renderer','painter'); clf
 
% Hidden states
%--------------------------------------------------------------------------
for p=1:2
    subplot(3,3,[1:3]+3*(p-1)),
    hax = gca;
    si    = spm_invNcdf(1 - 0.05);
    s     = [];
    if p == 1,
        xxfig = xxf;
        Sfig  = Sf;
        tit   = 'SCKF - forward pass';
    else
        xxfig = xxb;
        Sfig  = Sb;
        tit   = 'SCKS - backward pass';
    end
    for i = 1:T
        s = [s abs(diag(Sfig{i}))];
    end
 
    % conditional covariances
    %----------------------------------------------------------------------
    j           = [1:size(xxfig(:,:),1)];
    ss          = si*s(j,:);
    [ill indss] = sort(full(mean(ss,2)),'descend');
 
    pf = plot(1:T,xxfig,'linewidth',1.5);
    set(hax,'xlim',[1,T],'nextplot','add')
    box(hax,'on')
    for ic = 1:size(xxfig,1)
        col0 = get(pf(indss(ic)),'color');
        col  = (ones(1,3)-col0)*0.65 + col0;
        fill([(1:T) fliplr(1:T)],[(xxfig(indss(ic),:) + ss(indss(ic),:)) fliplr((xxfig(indss(ic),:) - ss(indss(ic),:)))],...
            'r',...
            'FaceColor',col,...
            'EdgeColor',col);
        hold on;
        COL{ic} = col0;
    end
    for ic = 1:size(xxfig,1)
        plot(xxfig(indss(ic),:),'color',COL{ic},'linewidth',0.75);
    end
    title(tit);
    grid(hax,'on');
    axis(hax,'tight');
end
 
 
% Parameter estimates
%--------------------------------------------------------------------------
subplot(3,3,7)
h = plot([1:length(ML)],ML);
if ~isempty(M(1).Q)
    AYlim = get(gca,'Ylim');
    if VBrun>=run
        bkg = ones(1,run)*max(AYlim);
    else
        bkg = [ones(1,VBrun)*max(AYlim),ones(1,abs(VBrun-run))*min(AYlim)];
    end
    a = area(bkg,min(AYlim));
    axis([1 run+1 AYlim(1) AYlim(2)]); hold on;
    set(a(1),'FaceColor',ones(1,3)*0.6,'EdgeColor',ones(1,3)*0.6)
    h = plot([1:length(ML)],ML);
    axis([1 run+1 AYlim(1) AYlim(2)]); hold on;
    if VBrun>=run
        text(run/2+1,mean(AYlim),'VB-SCKS','HorizontalAlignment','center','VerticalAlignment','top');
    else
        text(VBrun/2+1,mean(AYlim),'VB-SCKS','HorizontalAlignment','center','VerticalAlignment','top');
        text(VBrun+(abs(VBrun-run)/2)+1,mean(AYlim),'SCKS','HorizontalAlignment','center','VerticalAlignment','top');
    end
    set(gca,'Layer','top')
    hold off
end
set(h,'color','k','Marker','o','Markersize',4,'MarkerFaceColor','k','linewidth',1);
title('Log-Likelihood');
 
if ~isempty(ip)
    subplot(3,3,8)
    b1 = bar(ip,mean(xxb(wind,:),2)','FaceColor','k');
    set(b1,'BarWidth',0.5); hold on;
    title('Parameters');
end
 
% Standard deviation of residuals
%--------------------------------------------------------------------------
if ~isempty(M(1).Q)
    subplot(3,3,9)
    if VBrun >= run
        plot(RR'); hold on;
    else
        plot(RR'); hold on;
        switch(lower(M(1).Qf))
            case('all')
                sR  = RR';
            case('mean')
                sR  = repmat(mean(RR,2),1,T)';
            case('min')
                RRs = sort(RR,2,'descend');
                sR  = repmat(mean(RRs(:,round(T*0.90):end),2),1,T)';
            case('auto')
                dlim    = min(RR,[],2);
                ulim    = max(RR,[],2);
                if all(ulim./dlim<4)
                    sR  = RR';
                else
                    RRs  = sort(RR,2,'descend');
                    sR   = repmat(mean(RRs(:,round(T*0.90):end),2),1,T)';
                end
        end
        plot([1:T],sR,'r','linewidth',2);
    end
    axis(gca,'tight');
    title('Standard deviation of residuals');
    hold off;
end
drawnow

