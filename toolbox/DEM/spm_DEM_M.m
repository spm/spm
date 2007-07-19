function [M] = spm_dem_M(model,varargin)
% creates a [template] model structure
% FORMAT [M] = spm_dem_M(model,l,n)
% FORMAT [M] = spm_dem_M(model,X1,X2,...)
%
% model: "General linear model","GLM"
%        "Factor analysis","FA"
%        "Independent component analysis","ICA"
%        "Sparse coding","SC"
%        "State space model","SSM"
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
        if length(varargin);
            pE   = varargin;
        else
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

        end


        % non-hierarchical linear generative model (static)
        %==================================================================
    case{'general linear model','glm'}

        % Get design matrix
        %------------------------------------------------------------------
        if length(varargin);
            pE = varargin{1};
        else
            global SPM_DATA
            fig  = spm_data;
            set(fig,'Name','Please specify a design matrix')
            waitfor(fig)
            pE   = SPM_DATA;
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
            l    = varargin{1};
            l    = l(~~l);
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
            l          = varargin{1}(1);
        catch
            msgbox('please specify number of outputs')
            error(' ')
        end
        M          = spm_DEM_M('fa',[l l]);
        g          = length(M);

        % assume precisely small error at the first level
        %------------------------------------------------------------------
        M(1).hE    = 8;       
        M(1).hC    = 1/32;

        % assume unit causes
        %------------------------------------------------------------------
        M(2).V     = speye(M(2).l,M(2).l);

        % Independent component analysis (static - nonlinear)
        %==================================================================
    case{'independent component analysis','ica'}

        % create basic FA
        %------------------------------------------------------------------
        M          = spm_DEM_M('pca',varargin{:});
        g          = length(M);

        % and add supra-ordinate [nonlinear] level
        %------------------------------------------------------------------
        M(g + 1) = M(g);
        M(g).m   = M(g + 1).l;
        M(g).g   = 'v-tanh(v)*P';
        M(g).pE  = 0;
        M(g).pC  = 0;
        M(g).V   = speye(M(2).l,M(2).l)*1e6;


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


        % Nonlinear [Polynomial] factor analysis (static - nonlinear)
        %==================================================================
    case{'nonlinear model','nlfa'}

 

        % non-hierarchical linear generative model (dynamic)
        %===========================================================================
    case{'state space model','ssm'}

        % model specification
        %--------------------------------------------------------------------------
        f       = '[1; (1 + sin(P(2)*pi*x(1)) - P(1)*x(2) + exp(v))]';
        g       = '(x(2).^2)/5';
        M(1).x  = [1; 1];
        M(1).f  = inline(f,'x','v','P');
        M(1).g  = inline(g,'x','v','P');
        M(1).pE = [log(2) 0.04];
        M(1).V  = 1e4;

        M(2).v  = 0;
        M(2).V  = 2.4;

        % temporal correlations
        %--------------------------------------------------------------------------
        M(1).E.s  = 1/32;
        M(1).E.d  = 1;
        M(1).E.n  = 4;
        M(1).E.nD = 8;


    otherwise

        errordlg('unknown model; please add to spm_DEM_M')

end

% check and conplete model specification
%---------------------------------------------------------------------------
M       = spm_DEM_M_set(M);


