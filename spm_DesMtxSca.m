function varargout=spm_DesMtxSca(varargin)
% Scaling of design matrix portions to lie in [-1,1], for visual display.
% FORMAT [nX,nEnames]=spm_DesMtxSca(X1,X1names,X2,X2names,...);
%
% The arguments come in X,Enames pairs, where:
% X		- Design matrix, or a portion of one.
% Enames	- Names of effects represented by the columns of X, as
%                 given by spm_DesMtx
%
% Enames parameters can be omitted, the columns of the design matrix 
% portion are then normalised individually.
%
% nX		- Normalised image of the design matrix
% nEnames	- Names of the effects.
%_______________________________________________________________________
% 
% See also: spm_DesMtx
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%-Print warning of obsolescence
%-----------------------------------------------------------------------
warning('spm_DesMtxSca is grandfathered, use spm_DesMtx(''sca'',... instead')


%-Pass on arguments to spm_DesMtx
%-----------------------------------------------------------------------
varargout = cell(1,2);
[varargout{1:2}] = spm_DesMtx('sca',varargin{:});
