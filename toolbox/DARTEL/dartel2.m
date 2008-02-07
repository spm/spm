function varargout = dartel2(varargin)
% DARTEL 2D image registration stuff
%
% These functions are useful for testing out new ideas in 2D, but not so
% helpful for most practical applcations.
%_______________________________________________________________________
%
% FORMAT v = dartel2(v,g,f,param)
% v     - flow field n1*n2*2
% g     - image    n1*n2*n3
% f     - template n1*n2*n3
% param - 9 parameters (settings)
%         - [1] Regularisation type, can take values of
%           - 0 Linear elasticity
%           - 1 Membrane energy
%           - 2 Bending energy
%         - [2][3][4] Regularisation parameters
%           - For linear elasticity, the parameters
%             are mu, lambda, and id
%           - For membrane and bending energy, the parameters
%             are lambda, unused and id.
%         - [5] Levenberg-Marquardt regularisation
%         - [6] Number of Full Multigrid cycles
%         - [7] Number of relaxation iterations per cycle
%         - [8] K, such that 2^K time points are used to
%               generate the deformations.  A value of zero
%               indicates a small deformation model.
%         - [9] 0/1, indicating whether or not to use a
%               symmetric formulation for the objective function.
%
% This is for performing a single iteration of the DARTEL optimisation.
% Images can be scalar fields, in which case the objective function
% is the sum of squares difference.  Alternatively, images can be
% vector fields, in which case the objective function is the sum of squares
% difference between each scalar field + the sum of squares difference
% between one minus the sum of the scalar fields.
%
%_______________________________________________________________________
%
% FORMAT v = dartel2('cgs',A, b, param)
% v     - the solution
% A     - parameterisation of 2nd derivatives
% b     - parameterisation of first derivatives
% param - 6 parameters (settings)
%         - [1] Regularisation type, can take values of
%           - 0 Linear elasticity
%           - 1 Membrane energy
%           - 2 Bending energy
%         - [2][3] Pixel sizes
%         - [4][5][6] Regularisation parameters
%           - For Linear elasticity, the parameters
%             are mu, lambda, and id
%           - For membrane and bending energy, the parameters
%             are lambda, unused and id.
%         - [7] Tolerance.  Indicates required degree of accuracy.
%         - [8] Maximum number of iterations.
%
% This is for solving a set of equations using a conjugate gradient
% solver. This method is less efficient than the Full Multigrid.
% v = inv(A+H)*b
%
%_______________________________________________________________________
%
% FORMAT v = dartel2('fmg',A, b, param)
% v     - the solution n1*n2*2
% A     - parameterisation of 2nd derivatives 
% b     - parameterisation of first derivatives
% param - 6 parameters (settings)
%         - [1] Regularisation type, can take values of
%           - 0 Linear elasticity
%           - 1 Membrane energy
%           - 2 Bending energy
%         - [2][3] Pixel sizes
%         - [4][5][6] Regularisation parameters
%           - For linear elasticity, the parameters
%             are mu, lambda, and id
%           - For membrane and bending energy, the parameters
%             are lambda, unused and id.
%         - [7] Number of Full Multigrid cycles
%         - [8] Number of relaxation iterations per cycle
%
% Solve equations using a Full Multigrid method.  See Press et al
% for more information.
% v = inv(A+H)*b
%
%_______________________________________________________________________
%
% FORMAT y = dartel2('Exp', v)
% v - flow field.
% y - deformation field.
% A flow field is "exponentiated" to generate a deformation field
% using a scaling and squaring approach. See the work of Arsigny
% et al, or Moler's "19 dubious ways".
%
%_______________________________________________________________________
%
% FORMAT m = dartel2('vel2mom', v, param)
% v     - velocity (flow) field n1*n2*2.
% param - 4 parameters (settings)
%         - [1] Regularisation type, can take values of
%           - 0 Linear elasticity
%           - 1 Membrane energy
%           - 2 Bending energy
%         - [2][3] Pixel sizes
%         - [4][5][6] Regularisation parameters
%           - For linear elasticity, the parameters
%             are mu, lambda and id.
%           - For membrane and bending energy, the parameters
%             are lambda, unusaed and id.
% m       - `momentum' field n1*n2*2.
%
% Convert a flow field to a momentum field by m = H*v, where
% H is the large sparse matrix encoding some form of regularisation.
%
%_______________________________________________________________________
%
% FORMAT y3 = dartel2('comp',y1,y2)
% y1, y2 - deformation fields n1*n2*2.
% y3     - deformation field field n1*n2*2.
%
% Composition of two deformations y3 = y1(y2)
%
%_______________________________________________________________________
%
% FORMAT [y3,J3] = dartel3('comp', y1, y2, J1, J2)
% y1, y2 - deformation fields n1*n2*2.
% y3     - deformation field n1*n2*2.
% J1, J2 - Jacobian tensor fields n1*n2*2*2.
% J3     - Jacobian tensor field n1*n2*2*2.
%
% Composition of two deformations, with their Jacobian fields.
%
%_______________________________________________________________________
%
% FORMAT f2 = dartel2('samp', f1, y)
% f1 - input image(s) n1*n2*n3
% y  - points to sample n1*n2*2
% f2 - output image n1*n2*2
%
% Sample a function according to a deformation.
% f2 = f1(y)
%
%_______________________________________________________________________
%
% FORMAT v2 = dartel2('resize', v1, dim)
% v1  - input fields n1*n2*n3
% v2  - output field dim1*dim2*n3
% dim - output dimensions
%
% Resize a field according to dimensions dim.  This is
% a component of the FMG approach.
%
%_______________________________________________________________________
%
% FORMAT v3 = dartel2('brc', v1, v2)
% v1, v2, v3 - flow fields n1*n2*2
%
% Lie Bracket.  Useful for many things
% e.g. Baker-Campbell-Haussdorf series expansion.
% The Lie bracket is denoted by
% v3 = [v1,v2]
% and on scalar fields, is computed by
% v3 = J1*v2 - J2*v1, where J1 and J2 are the Jacobian
% tensor fields. For matrices, the Lie bracket is simply
% [A,B] = A*B-B*A
%
%_______________________________________________________________________
%
% Note that the boundary conditions are circulant throughout.
% Interpolation is trilinear, except for the resize function
% which uses a 2nd degree B-spline (without first deconvolving).
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: dartel2.m 1143 2008-02-07 19:33:33Z spm $

error('Not compiled for %s in MATLAB %s  (see make.m)\n', computer, version);

