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
% Ep        [Np x Nsamp] posterior density matrix. Parameter vector is of
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
% $Id: spm_dcm_bma.m 4067 2010-09-03 16:59:35Z maria $

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
if isfield(DCM,'a')
    dcm_fmri = 1;
    nreg = DCM.n;
    min  = DCM.M.m;
    dimD = 0;
else
    dcm_fmri = 0;
end

firstsub  = 1;
firstmod  = 1;

Ep  = [];
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

            params(i).model(kk).Ep  = subj(i).sess(1).model(sel).Ep;
            params(i).model(kk).vEp = spm_vec(params(i).model(kk).Ep);
            params(i).model(kk).Cp  = full(subj(i).sess(1).model(sel).Cp);

            if dcm_fmri
                dimDtmp = size(params(i).model(kk).Ep.D,3);
                if dimDtmp ~= 0, dimD = dimDtmp; firstsub = i; firstmod = kk;end
            end
            
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
                vEp         = Ep;
                Ep          = spm_unvec(Ep,pE);

                params(i).model(kk).Ep  = Ep;
                params(i).model(kk).vEp = vEp;
                params(i).model(kk).Cp  = Cp;

            end

            [evec, eval] = eig(params(i).model(kk).Cp);
            deig         = diag(eval);

            params(i).model(kk).dCp = deig;
            params(i).model(kk).vCp = evec;

        end   
    end
    
else % Use an FFX
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

            params(n).model(kk).Ep  = subj(n).sess(1).model(sel).Ep;
            params(n).model(kk).vEp = spm_vec(params(n).model(kk).Ep);
            params(n).model(kk).Cp  = full(subj(n).sess(1).model(sel).Cp);
            
            if dcm_fmri
                dimDtmp = size(params(n).model(kk).Ep.D,3);
                if dimDtmp ~= 0, dimD = dimDtmp; firstsub = n; firstmod = kk; end
            end

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
                vEp         = Ep;
                Ep          = spm_unvec(Ep,pE);

                params(n).model(kk).Ep  = Ep;
                params(n).model(kk).vEp = vEp;
                params(n).model(kk).Cp  = Cp;

            end

            [evec, eval] = eig(params(n).model(kk).Cp);
            deig         = diag(eval);

            params(n).model(kk).dCp = deig;
            params(n).model(kk).vCp = evec;
        end
    end

end

% Pre-allocate sample arrays
Np = length(params(firstsub).model(firstmod).vEp);

% get dimensions of a b c d parameters
if dcm_fmri

    Nr      = nreg*nreg;

    Etmp.A  = zeros(nreg,nreg,Nsamp);
    Etmp.B  = zeros(nreg,nreg,min,Nsamp);
    Etmp.C  = zeros(nreg,min,Nsamp);
    Etmp.D  = zeros(nreg,nreg,dimD,Nsamp);

    dima    = Nr;
    dimb    = Nr+Nr*min;
    dimc    = Nr+Nr*min+nreg*min;

end

clear Ep
disp('')
disp('Averaging models in Occams window...')

Ep_all = zeros(Np,Nsub);
Ep_sbj = zeros(Np,Nsub,Nsamp);
Ep     = zeros(Np,Nsamp);

for i=1:Nsamp
    % Pick a model
    if ~rfx
        m = spm_multrnd(post,1);
    end
    % Pick parameters from model for each subject
    for n=1:Nsub

        clear mu dsig vsig
        
        if rfx
            m = spm_multrnd(renorm(n).post,1);
        end

        mu                    = params(n).model(m).vEp;
        nmu                   = length(mu);
        dsig                  = params(n).model(m).dCp(1:nmu,1);
        vsig(:,:)             = params(n).model(m).vCp(1:nmu,1:nmu);

        tmp                   = spm_normrnd(mu,{dsig,vsig},1);
        
        Ep_all(1:nmu,n)      = tmp(:);
        Ep_sbj(1:nmu,n,i)    = Ep_all(1:nmu,n); 

    end
    
    % Average over subjects
    Ep(:,i) = mean(Ep_all,2);

end

% save mean parameters
Ep_avg     = mean(Ep,2);
Ep_std     = std(Ep,0,2);
Ep_avg     = spm_unvec(Ep_avg,params(1).model(1).Ep);
Ep_std     = spm_unvec(Ep_std,params(1).model(1).Ep);
bma.mEp    = Ep_avg;
bma.sEp    = Ep_std;

Ep_avgsbj  = mean(Ep_sbj,3);
Ep_stdsbj  = std(Ep_sbj,0,3);

for is=1:Nsub
    bma.mEps(is)=spm_unvec(Ep_avgsbj(:,is),params(1).model(1).Ep);
    bma.sEps(is)=spm_unvec(Ep_stdsbj(:,is),params(1).model(1).Ep);
end

if dcm_fmri
    bma.a  = spm_unvec(Ep(1:dima,:),Etmp.A);
    bma.b  = spm_unvec(Ep(dima+1:dimb,:),Etmp.B);
    bma.c  = spm_unvec(Ep(dimb+1:dimc,:),Etmp.C);
    if dimD ~=0
        bma.d  = spm_unvec(Ep(dimc+1:dimc+Nr*dimD,:),Etmp.D);
    else
        bma.d  = Etmp.D;
    end
end

% storing parameters
% -------------------------------------------------------------------------
bma.nsamp = Nsamp;
bma.oddsr = oddsr;
bma.Nocc  = Nocc;
bma.Mocc  = post_ind;
