function [BSD, BSDs, F, logBF] = spm_bsd_search(BSD, nograph)
%SPM_BSD_SUBMODELS Summary of this function goes here
%   Detailed explanation goes here
    try nograph; catch nograph = 0; end

    if isstruct(BSD)
        fqs = BSD.fqs; 
        Nm = length(fqs);
        fprintf('Evaluating %d possible models.\n', 2^Nm); 
    
        BSDs = cell(2^Nm,1);  
        F = zeros(2^Nm, 1); 
    
        BSDfull = spm_bsd(BSD); 
        BSD.options.separatenull = 0;
    
        parfor_progress(2^Nm);
    
        % parfor i = 1:2^Nm-1
        for i = 1:2^Nm-1

            tag = str2num(dec2bin((i-1) + 2^Nm)'); 
            tag = ~~tag(2:end); 
    
            BSDi = BSD; 
            BSDi.fqs = fqs(tag);
    
            BSDi.M.P.f = BSDfull.Ep.f(tag);
            BSDi.M.P.S = BSDfull.Ep.S(tag);
            BSDi.M.P.a = BSDfull.Ep.a(tag);
            BSDi.M.P.b = BSDfull.Ep.b;
    
            BSDi.name = spm_file(BSDi.name, 'suffix', sprintf('_search_%d', i)); 
            
            BSDi = spm_bsd(BSDi); 
            BSDs{i} = BSDi;
            F(i) = BSDi.F; 
    
            parfor_progress;
        end
    
        BSDs{end} = BSDfull;
        F(end) = BSDfull.F;
    else
        BSDs = BSD; 
        Nm = length(BSD); 
        F = cellfun(@(c) c.F, BSDs); 
        F = [F{:}]; 
    end


    Fm = repmat(F, 1, Nm);

    tag = arrayfun(@(a) str2num(a), dec2bin((1:2^Nm)-1));
    tag = ~~tag;

    [~, iF] = max(F); 
    BSD = BSDs{iF}; 

    logBF = zeros(Nm, 1); 
    for i = 1:Nm
        logBF(i) = mean(Fm(tag(:, i)), 1) - mean(Fm(~tag(:, i)), 1);
        % logBF(i) = logSumExp(Fm(tag(:,i)), 1) - logSumExp(Fm(~tag(:,i)), 1); 
    end


    if ~nograph 
        spm_figure('GetWin', 'BSD - Reduce'); 
        clf; 
        
        subplot(2,2,1);
        plot(F);
        title('free-energy')
        xlabel('models')
    
        subplot(2,2,2);
        bar(exp(F - logSumExp(F,1)));
        title('model posterior probability')
        xlabel('models');
    
        subplot(2,2,3);
    %     yregion(0, 20);
        bar(logBF);
        hold on; 
        yline(1, 'y:', 'LineWidth', 2);
        yline(3, 'k--', 'LineWidth',2);
        yline(5, 'r-', 'LineWidth', 2);
        hold off; 
        set(gca,'XTickLabel',fqs)
        ylim([0,20]);
        title('log Bayes factor')
        xlabel('frequencies'); 
    
        sel = logBF > 3; 
    
        P = zeros(sum(sel), 1); 
        [C, iA] = unique(tag(:, sel), 'rows');
        pos = sum(C, 2) > 0; 
        C = C(pos, :);
        iA = iA(pos, :); 
        [iA, iB] = sort(iA); 
        rtag = C(iB, :);
        Fm = Fm(iA, sel); 
        for i = 1:sum(sel)
            L = Fm(:, i);
            L(rtag(:, i)) = L(rtag(:, i));
            P(i) = exp(logSumExp(L(rtag(:,i)), 1) - logSumExp(L, 1)) ; 
        end
    
        
        subplot(2,2,4);
        % spm_plot_ci(BSD.Ep, BSD.fqs)
        bar(P);
        set(gca,'XTickLabel',fqs(sel))
        title('peak posterior probability')
        xlabel('frequencies'); 
    
%     parfor_progress(0);

    % spm_figure('GetWin', 'BSD - Reduce fit'); 
        spm_bsd_fit_summary(BSD); 
    end
end


function y = logSumExp(x, dim)
    X = max(x); 
    y = X + log(sum(exp(x-X), dim));
end
