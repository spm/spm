function bma = spm_dcm_bma(post,post_indx,subj,Nsamp,oddsr)
% Model-independent samples from DCM posterior  
% FORMAT bma = spm_dcm_bma (post,post_indx,subj,Nsamp,oddsr)
%
% post      [Ni x M] vector of posterior model probabilities
%           If Ni>1 then inference is based on subject-specific RFX posterior 
% post_indx models to use in BMA (position of models in subj sctructure)
% subj      subj(n).sess(s).model(m).fname: DCM filename
% Nsamp     Number of samples (default = 1e3)
% oddsr     posterior odds ratio for defining Occam's window (default=0, ie
%           all models used in average)
%
% theta     [Np x Nsamp] posterior density matrix. Parameter vector is of 
%           dimension Np and there are Nsamp samples
% Nocc      Number of models in Occam's window
%
%           For RFX BMA, different subject can have different models in
%           Occam's window (and different numbers of models in Occam's
%           window)
%
% This routine implements Bayesian averaging over models and subjects
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_dcm_bma.m 3734 2010-02-23 11:39:03Z maria $

if nargin < 4 || isempty(Nsamp)
    Nsamp = 1e3;
end
if nargin < 5 || isempty(oddsr)
    oddsr = 0;
end

Nsub = length(subj);
Nses = length(subj(1).sess);

% Number of regions
load(subj(1).sess(1).model(1).fname);
nreg = DCM.n;
min  = DCM.M.m;

theta  = [];
[Ni M] = size(post);

if Ni > 1 
    rfx = 1;
else
    rfx = 0;
end

if rfx
    
    for i = 1:Ni,
        
        mp          = max(post(i,:));
        post_ind{i} = find(post(i,:)>mp*oddsr);
        Nocc(i)     = length(post_ind{i});
        disp(' ');
        disp(sprintf('Subject %d has %d models in Occams window',i,Nocc(i)));
        
        if Nocc(i) == 0,
            return;
        end
        
        for occ = 1:Nocc(i),
            m = post_ind{i}(occ);
            disp(sprintf('Model %d, <p(m|Y>=%1.2f',m,post(i,m)));
        end
        
        % Renormalise post prob to Occam group
        renorm(i).post = post(i,post_ind{i});
        sp             = sum(renorm(i).post,2);
        renorm(i).post = renorm(i).post./(sp*ones(1,Nocc(i)));
        
        % Load DCM posteriors for models in Occam's window
        for kk = 1:Nocc(i),
            
            sel     = post_indx(post_ind{i}(kk));
            
            params(i).model(kk).Ep = subj(i).sess(1).model(sel).Ep;
            params(i).model(kk).Cp = full(subj(i).sess(1).model(sel).Cp);
              
            % Average sessions
            if Nses > 1
                
                clear miCp mEp
                disp('Averaging sessions...')
                
                for ss = 1:Nses 
                    
                    % Only parameters with non-zero prior variance
                    %------------------------------------------------------
                    sess_model.Cp = full(subj(i).sess(ss).model(sel).Cp);
                    pCdiag        = diag(full(sess_model.Cp));
                    wsel          = find(pCdiag);
                    
                    if ss == 1
                        wsel_first = wsel;
                    else
                        if ~(length(wsel) == length(wsel_first))
                            disp('Error: DCMs must have same structure');
                            return
                        end
                        if ~(wsel == wsel_first)
                            disp('Error: DCMs must have same structure');
                            return
                        end
                    end
 
                    % Get posterior precision matrix and mean
                    %------------------------------------------------------
                    Cp           = sess_model.Cp;
                    Ep           = spm_vec(subj(i).sess(ss).model(sel).Ep);
                    miCp(:,:,ss) = inv(full(Cp(wsel,wsel)));
                    mEp(:,ss)    = full(Ep(wsel));
 
                end
                
                % Average models using Bayesian fixed-effects analysis
                %==========================================================
                Cp(wsel,wsel) = inv(sum(miCp,3));
                
                pE            = subj(i).sess(ss).model(sel).Ep;
                weighted_Ep   = 0;
                for s = 1:Nses
                    weighted_Ep = weighted_Ep + miCp(:,:,s)*mEp(:,s);
                end
                Ep(wsel)    = Cp(wsel,wsel)*weighted_Ep;
                Ep          = spm_unvec(Ep,pE);
                
                params(i).model(kk).Ep = Ep;
                params(i).model(kk).Cp = Cp;
            
            end
            
            [evec, eval] = eig(params(i).model(kk).Cp);
            deig         = diag(eval);
            
            params(i).model(kk).dCp = deig;
            params(i).model(kk).vCp = evec;
            
        end
    end
else
    % Find models in Occam's window
    mp       = max(post);
    post_ind = find(post>mp*oddsr);
    Nocc     = length(post_ind);
    disp(' ');
    disp(sprintf('%d models in Occams window',Nocc));
    
    if Nocc == 0, return; end
    
    for occ = 1:Nocc,
        m = post_ind(occ);
        disp(sprintf('Model %d, p(m|Y)=%1.2f',m,post(m)));
    end
    
    % Renormalise post prob to Occam group
    post=post(post_ind);
    post=post/sum(post);
    
    % Load DCM posteriors for models in Occam's window
    for n=1:Nsub,
        
        for kk=1:Nocc,
            
            sel = post_indx(post_ind(kk));
            
            params(n).model(kk).Ep = subj(n).sess(1).model(sel).Ep;
            params(n).model(kk).Cp = full(subj(n).sess(1).model(sel).Cp);
            
            if Nses > 1
                
                clear miCp mEp
                disp('Averaging sessions...')
                
                % Average sessions
                for ss = 1:Nses
                    
                    % Only parameters with non-zero prior variance
                    %------------------------------------------------------
                    sess_model.Cp = full(subj(n).sess(ss).model(sel).Cp);
                    pCdiag        = diag(full(sess_model.Cp));
                    wsel          = find(pCdiag);
                    
                    if ss == 1
                        wsel_first = wsel;
                    else
                        if ~(length(wsel) == length(wsel_first))
                            disp('Error: DCMs must have same structure');
                            return
                        end
                        if ~(wsel == wsel_first)
                            disp('Error: DCMs must have same structure');
                            return
                        end
                    end
                    
                    % Get posterior precision matrix and mean
                    %------------------------------------------------------
                    Cp           = sess_model.Cp;
                    Ep           = spm_vec(subj(n).sess(ss).model(sel).Ep);
                    miCp(:,:,ss) = inv(full(Cp(wsel,wsel)));
                    mEp(:,ss)    = full(Ep(wsel));
                    
                end
                
                % Average models using Bayesian fixed-effects analysis
                %==========================================================
                Cp(wsel,wsel) = inv(sum(miCp,3));
                
                pE            = subj(n).sess(ss).model(sel).Ep;
                weighted_Ep   = 0;
                for s = 1:Nses
                    weighted_Ep = weighted_Ep + miCp(:,:,s)*mEp(:,s);
                end
                Ep(wsel)    = Cp(wsel,wsel)*weighted_Ep;
                Ep          = spm_unvec(Ep,pE);
                
                params(n).model(kk).Ep = Ep;
                params(n).model(kk).Cp = Cp;
                
            end
            
            [evec, eval] = eig(params(n).model(kk).Cp);
            deig         = diag(eval);
            
            params(n).model(kk).dCp = deig;
            params(n).model(kk).vCp = evec;
        end
    end
    
end

% Pre-allocate sample arrays
Nr        = nreg*nreg;
Np        = Nr + Nr*min + nreg*min + Nr*nreg;
theta_all = zeros(Np,Nsub);
nind      = 1;
mind      = 1;

disp('')
disp('Averaging models in Occams window...')
for i=1:Nsamp
    % Pick a model
    if ~rfx
        m = spm_multrnd(post,1);
    end
    % Pick parameters from model for each subject
    for n=1:Nsub
        
        if rfx
            m = spm_multrnd(renorm(n).post,1);
        end
        
        nD   = size(params(n).model(m).Ep.D,3);
        ndim = Nr + Nr*min + nreg*min + Nr*nD;
        
        if Np==ndim, nind=n; mind=m; end
        
        mu                  = zeros(Np,1);
        mu_tmp              = spm_vec(params(n).model(m).Ep);
        mu(1:ndim,1)        = mu_tmp(1:ndim,1);
       
        dsig                = zeros(Np,1);
        dsig(1:ndim,1)      = params(n).model(m).dCp(1:ndim,1);
        
        vsig                = zeros(Np,Np);
        vsig(1:ndim,1:ndim) = params(n).model(m).vCp(1:ndim,1:ndim);
        
        tmp                 = spm_normrnd(mu,{dsig,vsig},1);
        theta_all(:,n)      = tmp(:);
        
    end
    
    % Average over subjects
    if Nsub>1
        theta(:,i)       = mean(theta_all,2);
        theta_sbj(:,:,i) = theta_all;
    else
        theta(:,i)       = theta_all;
        theta_sbj(:,i)   = theta_all;
    end
    
    % unvec parameters
    % ---------------------------------------------------------------------
    Etmp.A = params(nind).model(mind).Ep.A;
    Etmp.B = params(nind).model(mind).Ep.B;
    Etmp.C = params(nind).model(mind).Ep.C;
    Etmp.D = params(nind).model(mind).Ep.D;
    
    Ep = spm_unvec(theta(:,i),Etmp);
    
    bma.a(:,:,i)    = Ep.A(:,:);
    bma.b(:,:,:,i)  = Ep.B(:,:,:);
    bma.c(:,:,i)    = Ep.C(:,:);
    bma.d(:,:,:,i)  = Ep.D(:,:,:);
    
end

% storing parameters
% -------------------------------------------------------------------------
bma.ma         = mean(bma.a,3);
bma.mb         = mean(bma.b,4);
bma.mc         = mean(bma.c,3);
bma.md         = mean(bma.d,4);

bma.theta      = theta;
bma.theta_sbj  = theta_sbj;
bma.mtheta_sbj = mean(theta_sbj,3);
bma.stheta_sbj = std(theta_sbj,0,3);

bma.nsamp      = Nsamp;
bma.oddsr      = oddsr;
bma.Nocc       = Nocc;
bma.Mocc       = post_ind;
