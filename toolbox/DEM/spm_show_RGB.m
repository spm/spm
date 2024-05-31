function [J,K] = spm_show_RGB(MDP,RGB,Nt,MOVIE)
% Graphical illustration of active inference
% FORMAT  [J,K] = spm_show_RGB(MDP,RGB,Nt,MOVIE)
% MDP   - Hierarchical (RG) generative model (inverted)
% RGB   - Structure for mapping discrete outcomes to RGB pixels
% Nt    - Number of timesteps to show images
% MOVIE - flag for movie [default: MOVIE = 1]
%
% [J,K] - predicted and observed RGB time-series
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
if nargin < 3, Nt    = 4; end                % Number of frames to show
if nargin < 4, MOVIE = 1; end

T     = size(O{1},2);
Nm    = numel(O);

% highest states and transitions
%--------------------------------------------------------------------------
subplot(Nm + 3,2,1)
image((1 - spm_cat(MDP.X))*64)
title(sprintf('Posterior (states) level %i',Nm),'FontSize',12)

% transition probabilities
%--------------------------------------------------------------------------
try
    B = spm_dir_norm(MDP.b{1});
catch
    B = MDP.B{1};
end

% image transition probabilities
%--------------------------------------------------------------------------
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
    title(sprintf('Discovered transitions (%i)',Nu),'FontSize',12)
    axis square
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
    % i = any(abs(diff(Q,1,2)) > 1/1024,2);
    % Q = Q(i,:);
    hold off, image((1 - Q)*64)
    title(sprintf('Predictive posterior (states) level %i',L),'FontSize',12)

    % uncomment for [relative] entropies
    %----------------------------------------------------------------------
    % R   = spm_dir_norm(R);
    % H   = sum(-R.*spm_log(R));
    % KL  = sum(R(:,2:end).*(spm_log(R(:,2:end)) - spm_log(R(:,1:end - 1))));
    % hold on, plot(KL,'r')
    % hold on, plot(10*H, 'g')

    % predicted paths (and entropy)
    %----------------------------------------------------------------------
    if L > 1
        subplot(Nm + 3,2,(n - 1)*2 + 4)
        Q = spm_cat(Y{L}(iE,:));
        % i = any(abs(diff(Q,1,2)) > 1/1024,2);
        % Q = Q(i,:);
        hold off, image((1 - Q)*64)
        title(sprintf('Predictive posterior (paths) level %i',L),'FontSize',12)

        % uncomment for [relative] entropies
        %----------------------------------------------------------------------
        % R   = spm_dir_norm(R);
        % H   = sum(-R.*spm_log(R));
        % KL  = sum(R(:,2:end).*(spm_log(R(:,2:end)) - spm_log(R(:,1:end - 1))));
        % hold on, plot(KL,'r')
        % hold on, plot(10*H, 'g')
    end

    % next level
    %----------------------------------------------------------------------
    if n < Nm
        mdp = mdp.MDP;
    end
end

% At the lowest (pixel) level
%--------------------------------------------------------------------------
I     = struct([]);
J     = uint8([]);
K     = uint8([]);
f     = 1;
for t = 1:T

    % return if figure is requested
    %----------------------------------------------------------------------
    if ~MOVIE && t > Nt
        return
    end

    % observed outcomes
    %----------------------------------------------------------------------
    X  = spm_O2rgb(O{1}(:,t),RGB);

    % predicted outcomes
    %----------------------------------------------------------------------
    P  = spm_O2rgb(Y{1}(:,t),RGB);

    % plot first level outcomes as a coloured image
    %----------------------------------------------------------------------
    subplot(Nm + 3,Nt,Nt*(Nm + 1) + min(Nt,t))
    spm_imshow(P)
    title(sprintf('Predicted: t = %i',t),'FontSize',10)


    % plot first level predictions as a coloured image
    %----------------------------------------------------------------------
    subplot(Nm + 3,Nt,Nt*(Nm + 2) + min(Nt,t))
    spm_imshow(X)
    if isfield(mdp,'S')
        title('Stimulus','FontSize',10)
    else
        title('Observed','FontSize',10)
    end

    % save movie
    %----------------------------------------------------------------------
    drawnow
    for i = 1:size(P,4)
        I(f).cdata    = P(:,:,:,i);
        I(f).colormap = [];
        J(:,:,:,f)    = P(:,:,:,i);
        K(:,:,:,f)    = X(:,:,:,i);
        f = f + 1;
    end
end

% Place movie in graphic subject
%--------------------------------------------------------------------------
try
    subplot(Nm + 3,Nt,Nt*(Nm + 2))
    h = get(gca,'Children');
    set(h(1),'Userdata',[])
    set(h(1),'Userdata',{I,8})
    set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')
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
