function [F,P] = spm_mvb_bmc
% multivariate Bayesian model comparison (Baysian decoding of a contrast)
% FORMAT [F,P] = spm_mvb_bmc
%__________________________________________________________________________


%-Get figure handles and set title
%--------------------------------------------------------------------------
Fmvb = spm_figure('GetWin','MVB');
if isempty(Fmvb)
    Fmvb = spm_figure('Create','MVB','Multivariate Bayes');
else
    clf
end

% get MVB results
%--------------------------------------------------------------------------
mvb  = spm_select(Inf,'mat','please select models',[],pwd,'MVB_*');
MVB  = load(deblank(mvb(1,:)));
MVB  = MVB.MVB;

% display
%==========================================================================
figure(Fmvb)

% if there is more than one MVB get maximum F for each
%--------------------------------------------------------------------------
if size(mvb,1) > 1
    
    X     = MVB.X;
    name  = {MVB.name(5:end)};
    F     = max(MVB.M.F);
    I     = 1;
    
    % check target is the same
    %----------------------------------------------------------------------
    for i = 2:size(mvb,1)
        MVB  = load(deblank(mvb(i,:)));
        MVB  = MVB.MVB;
        if ~any(X - MVB.X)
            name{end + 1} = MVB.name(5:end);
            F(end + 1)    = max(MVB.M.F);
            I(end + 1)    = i;
        end
    end

    % add null and model comparison
    %----------------------------------------------------------------------
    name{end + 1} = 'null';
    F(end + 1)    = MVB.M(1).F(1);
    F             = F - min(F);
    P             = exp(F - mean(F));
    P             = P./sum(P);
    
    % load best (non-null) model
    %----------------------------------------------------------------------
    [p i] = max(P(1:end - 1));
    MVB   = load(deblank(mvb(I(i),:)));
    MVB   = MVB.MVB;
    spm_mvb_display(MVB)
    
    % display model comparison
    %----------------------------------------------------------------------
    subplot(3,2,1)
    bar(F), hold on
    plot([0 length(F) + 1], [3 3],'r')
    plot([0 length(F) + 1],-[3 3],'r')
    plot([0 length(F) + 1], [5 5],'r:')
    plot([0 length(F) + 1],-[5 5],'r:')
    hold off
    set(gca,'XTickLabel',name);
    axis square
    grid on
    title({'log-evidence';'Model comparison'})

    subplot(3,2,2), cla
    text(.2,1/2,name','FontSize',12,'FontWeight','Bold')
    text(.6,1/2,num2str(P',' %.3f'),'FontSize',12,'FontWeight','Bold')
    axis square off
    title({'Posterior p-values';'Model comparison'})

else
    
    % display
    %----------------------------------------------------------------------
    spm_mvb_display(MVB)
    
    % probability relative to null
    %----------------------------------------------------------------------
    F     = MVB.M.F(2:end) - MVB.M.F(1);
    P     = exp(F);
    P     = P./(P + 1);
end


%-Reset title
%--------------------------------------------------------------------------
spm('Pointer','Arrow')



