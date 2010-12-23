function [M] = spm_dem_M(model,varargin)
% creates a [template] model structure
% FORMAT [M] = spm_dem_M(model,l,n)
% FORMAT [M] = spm_dem_M(model,X1,X2,...)
%
% model: 'General linear model','GLM'
%        'Factor analysis','FA'
%        'Independent component analysis','ICA'
%        'Sparse coding',"SC'
%        'convolution model'
%        'State space model','SSM',','Double Well'
%        'Lorenz'
%        'Ornstein_Uhlenbeck','OU'
%
% l(i) - number of outputs from level i
% n(i) - number of hidden states in level i
%
% Xi   - deisgn matrix for level i
%
%==========================================================================
% hierarchical generative model
%--------------------------------------------------------------------------
%   M(i).g  = y(t)  = g(x,v,P)    {inline function, string or m-file}
%   M(i).f  = dx/dt = f(x,v,P)    {inline function, string or m-file}
%
%   M(i).pE = prior expectation of p model-parameters
%   M(i).pC = prior covariances of p model-parameters
%   M(i).hE = prior expectation of h hyper-parameters (cause noise)
%   M(i).hC = prior covariances of h hyper-parameters (cause noise)
%   M(i).gE = prior expectation of g hyper-parameters (state noise)
%   M(i).gC = prior covariances of g hyper-parameters (state noise)
%   M(i).Q  = precision components (input noise)
%   M(i).R  = precision components (state noise)
%   M(i).V  = fixed precision (input noise)
%   M(i).W  = fixed precision (state noise)
%
%   M(i).m  = number of inputs v(i + 1);
%   M(i).n  = number of states x(i);
%   M(i).l  = number of output v(i);
%__________________________________________________________________________

switch lower(model)


    % hierarchical linear generative model - the basis for  static models
    % These subsume: Parametric emprical Bayes (PEB) models
    %                Mixed effects (MFX) models
    %======================================================================
    case{'hierarchical linear model','hlm','peb','mfx'}

        % Get design matrices for each level
        %------------------------------------------------------------------
        M(1).E.linear = 1;
        if ~isempty(varargin);
            pE   = varargin;
        else
            error('Please specify design matrices');
            
            global SPM_DATA;
            fig  = spm_data;
            set(fig,'Name','Please specify design matrices using ''next''')
            waitfor(fig)
            pE   = SPM_DATA;
        end

        for i = 1:length(pE)

            % level i
            %--------------------------------------------------------------
            M(i).m  = size(pE{i},2);
            M(i).l  = size(pE{i},1);
            M(i).g  = inline('P*v','x','v','P');
            M(i).pE = pE{i};
            M(i).hE = 0;
            M(i).hC = 16;

        end


        % non-hierarchical linear generative model (static)
        %==================================================================
    case{'general linear model','glm'}

        % Get design matrix
        %------------------------------------------------------------------
        if ~isempty(varargin);
            pE  = varargin{1};
        else
            global SPM_DATA
            fig = spm_data;
            set(fig,'Name','Please specify a design matrix')
            waitfor(fig)
            pE  = SPM_DATA;
        end

        % This is simply a one-level HLM with flat priors
        %------------------------------------------------------------------
        M       = spm_DEM_M('hlm',pE);


        % Factor analysis model (static) with unknown parameters
        %==================================================================
    case{'factor analysis','fa'}

        % Get orders (this routine will accommodate hierarchical FA)
        %------------------------------------------------------------------
        try
            l = varargin{1};
            l = l(~~l);
        catch
            errordlg('please specify number of inputs and ouputs')
        end

        % This is simply a HLM with unknown parameters
        %------------------------------------------------------------------
        g     = length(l) - 1;
        for i = 1:g
            pE{i} = sparse(l(i),l(i + 1));
            pC{i} = speye(l(i)*l(i + 1));
        end

        % create basic HLM and add prior uncertainty
        %------------------------------------------------------------------
        M     = spm_DEM_M('hlm',pE{:});
        for i = 1:g
            M(i).pC = pC{i};
        end


        % Principal component analysis (static - linear)
        % Equivalent to singular value decomposition (SVD)
        %==================================================================
    case{'principal component analysis','pca','svd'}

        % create basic FA
        %------------------------------------------------------------------
        try
            l    = varargin{1}(1);
        catch
            msgbox('please specify number of outputs')
            error(' ')
        end
        M        = spm_DEM_M('fa',[l l]);
        g        = length(M);

        % assume precisely small error at the first level
        %------------------------------------------------------------------
        M(1).hE  = 8;       
        M(1).hC  = 1/32;

        % assume unit causes
        %------------------------------------------------------------------
        M(2).V   = speye(M(2).l,M(2).l);

        % Independent component analysis (static - nonlinear)
        %==================================================================
    case{'independent component analysis','ica'}

        % create basic FA
        %------------------------------------------------------------------
        M        = spm_DEM_M('pca',varargin{:});
        g        = length(M);

        % and add supra-ordinate [nonlinear] level
        %------------------------------------------------------------------
        M(g + 1) = M(g);
        M(g).m   = M(g + 1).l;
        M(g).g   = 'v-tanh(v)*P';
        M(g).pE  = 0;
        M(g).pC  = 0;
        M(g).V   = speye(M(2).l,M(2).l)*exp(8);


        % Sparse coding, pICA (static - nonlinear)
        %==================================================================
    case{'sparse coding','pica','sc'}

        % create basic ica
        %------------------------------------------------------------------
        M          = spm_DEM_M('ica',varargin{:});
        g          = length(M);

        % and add a noise component to the lower levels
        %------------------------------------------------------------------
        for i = 1:(g - 2)
            M(i).hE = 1;
            M(i).V  = sparse(M(i).l,M(i).l);
        end

        % Linear convolution (State-space model)
        %==================================================================
    case{'convolution model'}
        
        % time-step
        %------------------------------------------------------------------
        try
            dt = varargin{1};
        catch
            dt = 1;
        end
        
        % smoothness
        %------------------------------------------------------------------
        M(1).E.linear = 1;                          % linear model
        M(1).E.s      = 1/2;                          % smoothness

        % level 1
        %------------------------------------------------------------------
        pE.f    = [-1  4;                           % the Jacobian for the
                   -2 -1]/(4*dt);                   % hidden sates
        pE.g    = [spm_dctmtx(4,2)]/4;              % the mixing parameters
        pE.h    = [1 0]';                           % input parameter
        M(1).n  = 2;
        M(1).f  = inline('P.f*x + P.h*v','x','v','P');
        M(1).g  = inline('P.g*x','x','v','P');
        M(1).pE = pE;                               % prior expectation

        % level 2
        %------------------------------------------------------------------
        M(2).l  = 1;                                % inputs
        M(2).V  = 1;                                % with shrinkage priors


        % Nonlinear [Polynomial] factor analysis (static - nonlinear)
        %==================================================================
    case{'nonlinear model','nlfa'}

 

        % non-hierarchical nonlinear generative model (dynamic)
        %==================================================================
    case{'state space model','ssm','double well'}
        
        
        % temporal correlations
        %------------------------------------------------------------------
        M(1).E.s  = 1/2;

        % model specification - 1st level
        %------------------------------------------------------------------
        f      = '(16*x./(1 + x.^2) - x/2 + v*2)/8';
        g      = '(x.^2)/16';
        M(1).x = 1;
        M(1).f = inline(f,'x','v','P');
        M(1).g = inline(g,'x','v','P');
        M(1).V = exp(2);
        M(1).W = exp(16);

        % model specification - 2nd level
        %------------------------------------------------------------------
        M(2).v = 0;
        M(2).V = 1/8;
        
                % non-hierarchical nonlinear generative model (dynamic)
        %==================================================================
    case{'lorenz'}
        
        
        % correlations
        %--------------------------------------------------------------------------
        M(1).E.linear = 3;
        M(1).E.s      = 1/8;

        % level 1
        %--------------------------------------------------------------------------
        P       = [18.0; -4.0; 46.92];
        x       = [0.9; 0.8; 30];
        f       = '[-P(1) P(1) 0; (P(3)-x(3)) -1 -x(1); x(2) x(1) P(2)]*x/32;';
        M(1).f  = f;
        M(1).g  = inline('sum(x)','x','v','P');
        M(1).x  = x;
        M(1).pE = P;
        M(1).V  = exp(0);
        M(1).W  = exp(16);

        % level 2
        %--------------------------------------------------------------------------
        M(2).v  = 0;
        M(2).V  = exp(16);
        

        % Ornstein_Uhlenbeck linear generative model (dynamic)
        %==================================================================
    case{'Ornstein_Uhlenbeck','ou'}
        
        % temporal correlations
        %------------------------------------------------------------------
        M(1).E.linear = 1;                          % linear model
        M(1).E.s  = 1/2;

        % model specification - 1st level
        %------------------------------------------------------------------
        f      = '-1/16*x + v';
        g      = 'x';
        M(1).x = 0;
        M(1).f = inline(f,'x','v','P');
        M(1).g = inline(g,'x','v','P');
        M(1).V = exp(8);
        M(1).W = exp(32);

        % model specification - 2nd level
        %------------------------------------------------------------------
        M(2).v = 0;
        M(2).V = exp(-8);
        

    otherwise

        errordlg('unknown model; please add to spm_DEM_M')
        return

end

% check and conplete model specification
%---------------------------------------------------------------------------
M       = spm_DEM_M_set(M);


