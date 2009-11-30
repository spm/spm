function [theta, Nocc] = spm_dcm_bma (post,post_indx,subj,Nsamp,oddsr)
% Model-independent samples from DCM posterior  
% FORMAT [theta, Nocc] = spm_dcm_bma (post,post_indx,subj,Nsamp,oddsr)
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
% $Id: spm_dcm_bma.m 3603 2009-11-30 18:56:50Z guillaume $

if nargin < 4 || isempty(Nsamp)
    Nsamp=1e3;
end
if nargin < 5 || isempty(oddsr)
    oddsr=0;
end

Nsub=length(subj);
Nses=length(subj(1).sess);

% Number of regions
load(subj(1).sess(1).model(1).fname);
nreg = DCM.n;
m    = DCM.M.m;
if isfield(DCM,'D')
    nonLin = 1;
else
    nonLin = 0;
end

theta=[];

[Ni M]=size(post);
if Ni > 1 
    rfx=1;
else
    rfx=0;
end

if rfx
    for i=1:Ni,
        mp=max(post(i,:));
        post_ind{i}=find(post(i,:)>mp*oddsr);
        Nocc(i)=length(post_ind{i});
        disp(' ');
        disp(sprintf('Subject %d has %d models in Occams window',i,Nocc(i)));
        if Nocc(i)==0,
            return;
        end
        
        for occ=1:Nocc(i),
            m=post_ind{i}(occ);
            disp(sprintf('Model %d, <p(m|Y>=%1.2f',m,post(i,m)));
        end
        
        % Renormalise post prob to Occam group
        renorm(i).post=post(i,post_ind{i});
        sp=sum(renorm(i).post,2);
        renorm(i).post=renorm(i).post./(sp*ones(1,Nocc(i)));
        
        % Load DCM posteriors for models in Occam's window
        for kk=1:Nocc(i),
            
            sel=post_indx(post_ind{i}(kk));
            
            params(i).model(kk).Ep=subj(i).sess(1).model(sel).Ep;
            params(i).model(kk).Cp=full(subj(i).sess(1).model(sel).Cp);
              
            % Average sessions
            if Nses > 1

                for ss=1:Nses 
                    clear miCp mEp

                    sess_model.Ep=subj(i).sess(ss).model(sel).Ep;
                    sess_model.Cp=full(subj(i).sess(ss).model(sel).Cp);
                                      
                    pCdiag = diag(sess_model.Cp);
                    wsel   = find(~(pCdiag==0));
                    
                    if nonLin
                       % but keep the D values if present !
                       npABC = n*n + n*n*m + n*m + 1 ; % nr of parameters in A,B,C+1
                       cwsel = wsel; cwsel(max(find(wsel<=npABC))+(1:6*n))=[];
                    else
                    cwsel = wsel(1:end-6*nreg);
                    end

                    % Get posterior precision matrix from model
                    miCp(:,:,ss) = inv(full(sess_model.Cp(cwsel,cwsel)));
                    % Get posterior mean from model
                    mEp(:,ss)    = full(sess_model.Ep(cwsel));

                end
                
                % Average models using Bayesian fixed effects analysis -> average Ep,Cp
                % averaged posterior covariance
                final_iCp = sum(miCp,3);
                Cp        = inv(final_iCp);
                % averaged posterior mean
                weighted_Ep = zeros(length(cwsel),1);
                for ss = 1:Nses,
                    weighted_Ep = weighted_Ep + miCp(:,:,ss)*mEp(:,ss);
                end
                Ep = Cp*weighted_Ep;

                params(i).model(kk).Ep(cwsel)=Ep;
                params(i).model(kk).Cp(cwsel,cwsel)=Cp;
                
            end
            
            [evec, eval] = eig(params(i).model(kk).Cp);
            deig=diag(eval);
            
            params(i).model(kk).dCp=deig;
            params(i).model(kk).vCp=evec;
            
        end
    end
else
    % Find models in Occam's window
    mp=max(post);
    post_ind=find(post>mp*oddsr);
    Nocc=length(post_ind);
    disp(' ');
    disp(sprintf('%d models in Occams window',Nocc));
    if Nocc==0, return; end
    for occ=1:Nocc,
        m=post_ind(occ);
        disp(sprintf('Model %d, p(m|Y)=%1.2f',m,post(m)));
    end
    
    % Renormalise post prob to Occam group
    post=post(post_ind);
    post=post/sum(post);
    
    % Load DCM posteriors for models in Occam's window
    for n=1:Nsub,
        for kk=1:Nocc,
            sel=post_indx(post_ind(kk));
           
            params(n).model(kk).Ep=subj(n).sess(1).model(sel).Ep;
            params(n).model(kk).Cp=full(subj(n).sess(1).model(sel).Cp);
              
            % Average sessions
            if Nses > 1
                
                for ss=1:Nses 
                    clear miCp mEp
                    
                    sess_model.Ep=subj(n).sess(ss).model(sel).Ep;
                    sess_model.Cp=full(subj(n).sess(ss).model(sel).Cp);
                                      
                    pCdiag = diag(sess_model.Cp);
                    wsel   = find(~(pCdiag==0));
                    
                    if nonLin
                       % but keep the D values if present !
                       npABC = n*n + n*n*m + n*m + 1 ; % nr of parameters in A,B,C+1
                       cwsel = wsel; cwsel(max(find(wsel<=npABC))+(1:6*n))=[];
                    else
                    cwsel = wsel(1:end-6*nreg);
                    end

                    % Get posterior precision matrix from model
                    miCp(:,:,ss) = inv(full(sess_model.Cp(cwsel,cwsel)));
                    % Get posterior mean from model
                    mEp(:,ss)    = full(sess_model.Ep(cwsel));

                end
                
                % Average models using Bayesian fixed effects analysis -> average Ep,Cp
                % averaged posterior covariance
                final_iCp = sum(miCp,3);
                Cp        = inv(final_iCp);
                % averaged posterior mean
                weighted_Ep = zeros(length(cwsel),1);
                for ss = 1:Nses,
                    weighted_Ep = weighted_Ep + miCp(:,:,ss)*mEp(:,ss);
                end
                Ep = Cp*weighted_Ep;

                params(n).model(kk).Ep(cwsel)=Ep;
                params(n).model(kk).Cp(cwsel,cwsel)=Cp;
                
            end 
            
            [evec, eval] = eig(params(n).model(kk).Cp);
            deig=diag(eval);
            
            params(n).model(kk).dCp=deig;
            params(n).model(kk).vCp=evec;
        end
    end
end

% Pre-allocate sample arrays
Np = length(spm_vec(params(1).model(1).Ep));
theta_all = zeros(Np,Nsub);

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
        mu   = params(n).model(m).Ep;
        mu   = spm_vec(mu);
        sig  = params(n).model(m).Cp;
        dsig = params(n).model(m).dCp;
        vsig = params(n).model(m).vCp;
        tmp  = spm_samp_gauss(mu,{dsig,vsig},1);
        theta_all(:,n) = tmp(:);
    end
    
    % Average over subjects
    if Nsub>1
        theta(:,i) = mean(theta_all,2);
    else
        theta(:,i) = theta_all;
    end
end

