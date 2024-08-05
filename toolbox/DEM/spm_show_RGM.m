function spm_show_RGM(MDP)
% Graphical illustration of active inference
% FORMAT spm_show_RGB(MDP,RGB,Nt,MOVIE)
% MDP   - Hierarchical (RG) generative model (inverted)
% Nt    - Number of timesteps to show images
% MOVIE - flag for movie [default: MOVIE = 1]
%
%--------------------------------------------------------------------------
% This auxiliary routine plots the hierarchal outcomes from a deep
% generative model. The separation of temporal scales is illustrated by
% showing the outcomes as a movie saved as the user data in the current
% figure. This routine separates predicted states and paths and shows the
% first few (predictions and observations) in pixel space.
%__________________________________________________________________________

% Get observations and predictions
%--------------------------------------------------------------------------
O     = spm_get_O(MDP);                      % observations
try
    Y = MDP.Q.Y;
catch
    Y = {};
end
Y{end + 1} = MDP.Y;                          % predictions

% Number of timesteps and levels
%==========================================================================
if nargin < 2, Nt    = 4; end                % Number of frames to show
if nargin < 3, MOVIE = 1; end

T     = size(O{1},2);
Nm    = numel(O);

% highest states and transitions
%--------------------------------------------------------------------------
subplot(Nm + 3,2,1)
if size(MDP.X,1) == 1
    image((1 - spm_cat(MDP.X(:)))*64)
else
    image((1 - spm_cat(MDP.X))*64)
end
title(sprintf('Posterior (states) level %i',Nm),'FontSize',12)

% transition probabilities
%--------------------------------------------------------------------------
Nf   = spm_MDP_size(MDP);
if Nf > 1

    % for each factor
    %----------------------------------------------------------------------
    for f = 1:Nf
        try
            B = spm_dir_norm(MDP.b{f});
        catch
            B = MDP.B{f};
        end
        subplot(Nm + 3,2*Nf,Nf + f)
        B = sum(B,3) > 1/16;
        if size(B,1) > 128
            spm_spy(B,8);
        else
            imagesc(1 - B)
        end
        title(sprintf('Transitions (f = %i)',f),'FontSize',12)
        axis square
    end

else

    % image transition probabilities
    %----------------------------------------------------------------------
    try
        B = spm_dir_norm(MDP.b{1});
    catch
        B = MDP.B{1};
    end
    Nu    = size(B,3);
    if Nu < 4
        for u = 1:Nu
            subplot(Nm + 3,2*Nu,Nu + u)
            if size(B,1) > 128
                spm_spy(B(:,:,u) > 1/16,8);
            else
                imagesc(1 - B(:,:,u))
            end
            title(sprintf('Transitions (u = %i)',u),'FontSize',12)
            axis square
        end
    else
        subplot(Nm + 3,2,2)
        B = sum(B,3) > 1/16;
        if size(B,1) > 128
            spm_spy(B,8);
        else
            imagesc(1 - B)
        end
        title(sprintf('Transitions (Nu = %i)',Nu),'FontSize',12)
        axis square
    end
end

% preditive posteriors over hidden states and paths: raster format
%--------------------------------------------------------------------------
mdp = MDP;
for n = 1:Nm
    
    L   = Nm - n + 1;                        % hierarchical level
    if L > 1
        iD  = spm_cat(mdp.MDP.id.D);         % states
        iE  = spm_cat(mdp.MDP.id.E);         % paths
    else
        iD  = 1:size(Y{L},1);
        iE  = [];
    end


    % predicted states (and entropy)
    %----------------------------------------------------------------------
    subplot(Nm + 3,2,(n - 1)*2 + 3)
    Q = spm_cat(Y{L}(iD,:));
    hold off, image((1 - Q)*64)
    title(sprintf('Predictive posterior (states) level %i',L),'FontSize',12)

    % predicted paths (and entropy)
    %----------------------------------------------------------------------
    if L > 1
        subplot(Nm + 3,2,(n - 1)*2 + 4)
        Q = spm_cat(Y{L}(iE,:));
        hold off, image((1 - Q)*64)
        title(sprintf('Predictive posterior (paths) level %i',L),'FontSize',12)

    end

    % next level
    %----------------------------------------------------------------------
    if n < Nm
        mdp = mdp.MDP;
    end
end

% At the lowest level
%==========================================================================

% observed outcomes
%--------------------------------------------------------------------------
X  = (1 - spm_cat(O{1}));

% predicted outcomes
%--------------------------------------------------------------------------
P  = (1 - spm_cat(Y{1}));

% first level outcomes
%--------------------------------------------------------------------------
subplot(Nm + 3,1,Nm + 2)
imagesc(P)
title('Predicted','FontSize',10)

% first level predictions
%-------------------------------------------------------------------------
subplot(Nm + 3,1,Nm + 3)
imagesc(X)
if isfield(mdp,'S')
    title('Stimulus','FontSize',10)
else
    title('Observed','FontSize',10)
end

return


% subroutines
%==========================================================================

function Q = spm_get_O(RDP)
% hierarchical outcomes from a recursive MDP
% FORMAT Q = spm_get_O(RDP)
%
% Extracts hierarchal outcomes from a recursive MDP that have been
% accumulated during recursive inversion
%__________________________________________________________________________

% recover outcomes
%--------------------------------------------------------------------------
try
    Q = RDP.Q.O;
catch
    Q = {};
end
Q{end + 1} = RDP.O;

return
