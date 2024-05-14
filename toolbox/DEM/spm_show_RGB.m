function spm_show_RGB(MDP,RGB,Nt,MOVIE)
% Graphical illustration of active inference
% FORMAT spm_show_MDP(MDP,RGB,[Nt],[FIG])
% MDP   - Hierarchical (RG) generative model (inverted)
% RGB   - Structure for mapping discrete outcomes to RGB pixels
% Number of timesteps to show images
% MOVIE - flag for figure [default: MOVIE = 1]
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

subplot(Nm + 3,2,2)
imagesc(1 - MDP.B{1}), axis square
title('Transitions','FontSize',12)


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
    hold off, image((1 - spm_cat(Y{L}(iD,:)))*64)
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
        hold off, image((1 - spm_cat(Y{L}(iE,:)))*64)
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
    title(sprintf('Predictons: t = %i',t),'FontSize',10)


    % plot first level predictions as a coloured image
    %----------------------------------------------------------------------
    if isfield(mdp,'S')
        subplot(Nm + 3,Nt,Nt*(Nm + 2) + min(Nt,t))
        spm_imshow(X)
        title(sprintf('Stimulus: t = %i',t),'FontSize',10)
    end

    % save movie
    %----------------------------------------------------------------------
    drawnow
    for i = 1:size(P,4)
        I(f).cdata = P(:,:,:,i);
        I(f).colormap  = [];
        f = f + 1;
    end

end

% Place movie in graphic subject
%--------------------------------------------------------------------------
h = get(gca,'Children');
set(h(1),'Userdata',[])
set(h(1),'Userdata',{I,8})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')

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
