function b = spm_fzero(varargin)
% Find a zero of a function of one variable
% FORMAT b = spm_fzero(FunFcn,x,tol,trace,varargin)
% FunFcn - String containing name of function
% x      - initial guess
% tol    - (optional) relative tolerance for the convergence test
%          (leave empty for default value)
% trace  - (optional) non-zero value triggers a printing trace of the steps
%          (leave empty for default value)
% p1-p9  - additional constant parameters, passed to function FunFcn.
%          spm_fzero finds zero of eval([FunFcn '(x,p1,p2,..)']) over
%          values of x.
%__________________________________________________________________________
%
% spm_fzero is obsolete in Matlab5, its functionality is superceeded by
% the fzero command in the core Matlab FunFun library.
%
% spm_fzero passes the arguments on to fzero, warning of its obsolescence.
%
% ( In MatLab4, fzero did not allow additional parameters to be passed   )
% ( to the target function, and spm_fzero was created from the fzero     )
% ( source code to allow this.                                           )
%__________________________________________________________________________
% %W% Andrew Holmes %E%

%-Print warning of obsolescence
%-----------------------------------------------------------------------
warning('spm_fzero is grandfathered in ML5, using fzero instead')


%-Pass on arguments to MatLab5 fzero
%-----------------------------------------------------------------------
b = fzero(varargin{:});
