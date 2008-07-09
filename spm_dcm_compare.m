function spm_dcm_compare(P)
% Compare two or more estimated models
% FORMAT spm_dcm_compare(P)
%
% P  - a char or cell array of DCM filenames
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Klaas Enno Stephan
% $Id: spm_dcm_compare.m 1900 2008-07-09 14:57:34Z guillaume $


% Get DCM filenames
%--------------------------------------------------------------------------
if ~nargin
    P = spm_select([2 inf],'^DCM.*\.mat$','select DCM*.mat files');
end
if ~iscell(P), P = cellstr(P); end
num_models = numel(P);

% Load anc check all models and compute their evidence
%--------------------------------------------------------------------------
for model_index=1:num_models
    
    try
        load(P{model_index},'DCM','-mat');
    catch
        error('Cannot load DCM file "%s".',P{model_index});
    end

    % Check that none of the models is an averaged DCM
    %----------------------------------------------------------------------
    if isfield(DCM,'averaged') && DCM.averaged
        str = {...
            ['Model ' P{model_index} ' is an averaged DCM. '],...
             'Please note that model comparison is not valid for averaged DCMs. ',...
             'Procedure aborted.'};
        spm('alert*',str,'DCM <Compare>',spm('Cmdline'));
        return
    end
    
    % Check that all models refer to the same set of VOIs
    %----------------------------------------------------------------------
    if model_index == 1
        VOIs = {DCM.xY.name};
    else
        if ~isequal(VOIs,{DCM.xY.name})
            str = {...
                'Selected models contain different sets of VOIs!',...
                'Please note that model comparison is only valid for models with identical time series (VOIs).',...
                'Procedure aborted.'};
            spm('alert*',str,'DCM <Compare>',spm('Cmdline'));
            return
        end
    end
    
    % Compute Model Evidence
    %----------------------------------------------------------------------
    evidence(model_index) = spm_dcm_evidence(DCM);
    aic(model_index)      = evidence(model_index).aic_overall;
    bic(model_index)      = evidence(model_index).bic_overall;
    
end


% Get and plot posterior probabilities of models assuming each model has 
% equal probability a-priori
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
figure(Fgraph);
spm_clf(Fgraph);

% Posterior probabilities of models from AIC
%--------------------------------------------------------------------------
maic = aic - mean(aic);
pm   = exp(maic) / sum(exp(maic));

subplot(2,1,1);
bar(pm);
Vscale = axis;
Vscale(4) = 1;
axis(Vscale);
set(gca,'FontSize',18);
ylabel('p(y|m)');  % given the flat prior p(m), p(y|m)=p(m|y)
xlabel('m');
title('Posterior probabilities of models from AIC');

% Posterior probabilities of models from BIC
%--------------------------------------------------------------------------
mbic = bic - mean(bic);
pm   = exp(mbic) / sum(exp(mbic));

subplot(2,1,2);
bar(pm);
Vscale = axis;
Vscale(4) = 1;
axis(Vscale);
set(gca,'FontSize',18);
ylabel('p(y|m)');  % given the flat prior p(m), p(y|m)=p(m|y)
xlabel('m');
title('Posterior probabilities of models from BIC');

% Output model comparison details to MATLAB command window
%--------------------------------------------------------------------------
for ii=1:num_models
    for jj=1:num_models
        if jj ~= ii

            disp(repmat('-',1,72));
            disp(sprintf('Model %d: %s',ii,P{ii}));
            disp('          versus ');
            disp(sprintf('Model %d: %s',jj,P{jj}));
            diff = 1;
            disp(' ');
            disp('All costs are in units of binary bits.');
            disp(' ');
            for k=1:size(DCM.A,1),
                nats = -(evidence(ii).region_cost(k) - evidence(jj).region_cost(k));
                bits = nats / log(2);
                disp(sprintf('Region %s: relative cost  = %1.4g, BF= %1.4g',DCM.Y.name{k},bits,2^(-bits)));
            end
            % AIC penalty
            %--------------------------------------------------------------
            nats = evidence(ii).aic_penalty - evidence(jj).aic_penalty;
            bits = nats / log(2);
            disp(sprintf('AIC Penalty = %1.4f, BF = %1.4g',bits,2^(-bits)));

            % BIC penalty
            %--------------------------------------------------------------
            nats = evidence(ii).bic_penalty - evidence(jj).bic_penalty;
            bits = nats / log(2);
            disp(sprintf('BIC Penalty = %1.4f, BF = %1.4g',bits,2^(-bits)));

            % AIC overall
            %--------------------------------------------------------------
            nats = -diff * (aic(ii)-aic(jj));
            bits = nats / log(2);
            bf_aic = 2^(-bits);
            disp(sprintf('AIC Overall = %1.4f, BF = %1.4g',bits,bf_aic));

            % BIC overall
            %--------------------------------------------------------------
            nats = -diff * (bic(ii)-bic(jj));
            bits = nats / log(2);
            bf_bic = 2^(-bits);
            disp(sprintf('BIC Overall = %1.4f, BF = %1.4g',bits,bf_bic));
            disp(' ');

            % Evaluate results
            %--------------------------------------------------------------
            if ((bf_bic > 1) && (bf_aic < 1)) || ((bf_bic < 1) && (bf_aic > 1))
                % AIC and BIC do not concur
                disp('AIC and BIC disagree about which model is superior - no decision can be made.');
            else
                % AIC and BIC do concur - assign statement according to Raftery classification
                raftery_labels = {'Weak','Positive','Strong','Very strong'};
                raftery_thresh = [  3         20      150         inf];
                if bf_aic > 1
                    minBF          = min(bf_aic,bf_bic);
                    posBF          = find(sort([raftery_thresh minBF])==minBF);
                    m              = ii;
                else
                    raftery_thresh = fliplr(1./raftery_thresh);
                    minBF          = max(bf_aic,bf_bic);
                    posBF          = find(sort(-1*[raftery_thresh minBF])==-minBF);
                    m              = jj;
                end
                label          = raftery_labels{posBF};
                disp(sprintf([label ' evidence in favour of model %d'],m));
                disp(sprintf('Bayes factor >= %1.4g', minBF));
                disp(' ');
                disp(repmat('-',1,72));
            end

        end
    end
end

spm_input('Model comparison details in MATLAB command window',1,'d');