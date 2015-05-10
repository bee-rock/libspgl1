function [x,r,g,info] = spgl1( A, b, tau, sigma, x, options )

%SPGL1  Solve basis pursuit, basis pursuit denoise, and LASSO
%
% [x, r, g, info] = spgl1(A, b, tau, sigma, x0, options)
%
% ---------------------------------------------------------------------
% Solve the basis pursuit denoise (BPDN) problem
%
% (BPDN)   minimize  ||x||_1  subj to  ||Ax-b||_2 <= sigma,
%
% or the l1-regularized least-squares problem
%
% (LASSO)  minimize  ||Ax-b||_2  subj to  ||x||_1 <= tau.
% ---------------------------------------------------------------------
%
% INPUTS
% ======
% A        is an m-by-n matrix, explicit or an operator.
%          If A is a function, then it must have the signature
%
%          y = A(x,mode)   if mode == 1 then y = A x  (y is m-by-1);
%                          if mode == 2 then y = A'x  (y is n-by-1).
%
% b        is an m-vector.
% tau      is a nonnegative scalar; see (LASSO).
% sigma    if sigma != inf or != [], then spgl1 will launch into a
%          root-finding mode to find the tau above that solves (BPDN).
%          In this case, it's STRONGLY recommended that tau = 0.
% x0       is an n-vector estimate of the solution (possibly all
%          zeros). If x0 = [], then SPGL1 determines the length n via
%          n = length( A'b ) and sets  x0 = zeros(n,1).
% options  is a structure of options from spgSetParms. Any unset options
%          are set to their default value; set options=[] to use all
%          default values.
%
% OUTPUTS
% =======
% x        is a solution of the problem
% r        is the residual, r = b - Ax
% g        is the gradient, g = -A'r
% info     is a structure with the following information:
%          .tau     final value of tau (see sigma above)
%          .rNorm   two-norm of the optimal residual
%          .rGap    relative duality gap (an optimality measure)
%          .gNorm   Lagrange multiplier of (LASSO)
%          .stat    = 1 found a BPDN solution
%                   = 2 found a BP sol'n; exit based on small gradient
%                   = 3 found a BP sol'n; exit based on small residual
%                   = 4 found a LASSO solution
%                   = 5 error: too many iterations
%                   = 6 error: linesearch failed
%                   = 7 error: found suboptimal BP solution
%                   = 8 error: too many matrix-vector products
%          .time    total solution time (seconds)
%          .nProdA  number of multiplications with A
%          .nProdAt number of multiplications with A'
%
% OPTIONS
% =======
% Use the options structure to control various aspects of the algorithm:
%
% options.fid         File ID to direct log output
%        .verbosity   0=quiet, 1=some output, 2=more output.
%        .iterations  Max. number of iterations (default if 10*m).
%        .bpTol       Tolerance for identifying a basis pursuit solution.
%        .optTol      Optimality tolerance (default is 1e-4).
%        .decTol      Larger decTol means more frequent Newton updates.
%        .subspaceMin 0=no subspace minimization, 1=subspace minimization.
%
% EXAMPLE
% =======
%   m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
%   p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
%   A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
%   b  = A*x0 + 0.005 * randn(m,1);
%   opts = spgSetParms('optTol',1e-4);
%   [x,r,g,info] = spgl1(A, b, 0, 1e-3, [], opts); % Find BP sol'n.
%
% AUTHORS
% =======
%  Ewout van den Berg (ewout78@cs.ubc.ca)
%  Michael P. Friedlander (mpf@cs.ubc.ca)
%    Scientific Computing Laboratory (SCL)
%    University of British Columbia, Canada.
%
% BUGS
% ====
% Please send bug reports or comments to
%            Michael P. Friedlander (mpf@cs.ubc.ca)
%            Ewout van den Berg (ewout78@cs.ubc.ca)

% 15 Apr 07: First version derived from spg.m.
%            Michael P. Friedlander (mpf@cs.ubc.ca).
%            Ewout van den Berg (ewout78@cs.ubc.ca).
% 17 Apr 07: Added root-finding code.
% 18 Apr 07: sigma was being compared to 1/2 r'r, rather than
%            norm(r), as advertised.  Now immediately change sigma to
%            (1/2)sigma^2, and changed log output accordingly.
% 24 Apr 07: Added quadratic root-finding code as an option.
% 24 Apr 07: Exit conditions need to guard against small ||r||
%            (ie, a BP solution).  Added test1,test2,test3 below.
% 15 May 07: Trigger to update tau is now based on relative difference
%            in objective between consecutive iterations.
% 15 Jul 07: Added code to allow a limited number of line-search
%            errors. 
% 23 Feb 08: Fixed bug in one-norm projection using weights. Thanks
%            to Xiangrui Meng for reporting this bug.
% 26 May 08: The simple call spgl1(A,b) now solves (BPDN) with sigma=0.
% 18 Mar 13: Reset f = fOld if curvilinear line-search fails.
%            Avoid computing the Barzilai-Borwein scaling parameter
%            when both line-search algorithms failed.
% 09 Sep 13: Recompute function information at new x when tau decreases.
%            Fixed bug in subspace minimization. Thanks to Yang Lei
%            for reporting this bug.
    
%   ----------------------------------------------------------------------
%   This file is part of SPGL1 (Spectral Projected-Gradient for L1).
%
%   Copyright (C) 2007 Ewout van den Berg and Michael P. Friedlander,
%   Department of Computer Science, University of British Columbia, Canada.
%   All rights reserved. E-mail: <{ewout78,mpf}@cs.ubc.ca>.
%
%   SPGL1 is free software; you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as
%   published by the Free Software Foundation; either version 2.1 of the
%   License, or (at your option) any later version.
%
%   SPGL1 is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
%   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General
%   Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with SPGL1; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
%   USA
%   ----------------------------------------------------------------------
REVISION = '$Revision: 1017 $';
DATE     = '$Date: 2008-06-16 22:43:07 -0700 (Mon, 16 Jun 2008) $';
REVISION = REVISION(11:end-1);
DATE     = DATE(35:50);

tic;              % Start your watches!
m = length(b);

%----------------------------------------------------------------------
% Check arguments. 
%----------------------------------------------------------------------
if ~exist('options','var'), options = []; end
if ~exist('x','var'), x = []; end
if ~exist('sigma','var'), sigma = []; end
if ~exist('tau','var'), tau = []; end

if nargin < 2 || isempty(b) || isempty(A)
   error('At least two arguments are required');
elseif isempty(tau) && isempty(sigma)
   tau = 0;
   sigma = 0;
   singleTau = false;
elseif isempty(sigma) % && ~isempty(tau)  <-- implied
   singleTau = true;
else
   if isempty(tau)
      tau = 0;
   end
   singleTau = false;
end

%----------------------------------------------------------------------
% Grab input options and set defaults where needed. 
%----------------------------------------------------------------------
defaultopts = spgSetParms(...
'fid'        ,      1 , ... % File ID for output
'verbosity'  ,      2 , ... % Verbosity level
'iterations' ,   10*m , ... % Max number of iterations
'nPrevVals'  ,      3 , ... % Number previous func values for linesearch
'bpTol'      ,  1e-06 , ... % Tolerance for basis pursuit solution 
'lsTol'      ,  1e-06 , ... % Least-squares optimality tolerance
'optTol'     ,  1e-04 , ... % Optimality tolerance
'decTol'     ,  1e-04 , ... % Req'd rel. change in primal obj. for Newton
'stepMin'    ,  1e-16 , ... % Minimum spectral step
'stepMax'    ,  1e+05 , ... % Maximum spectral step
'rootMethod' ,      2 , ... % Root finding method: 2=quad,1=linear (not used).
'activeSetIt',    Inf , ... % Exit with EXIT_ACTIVE_SET if nnz same for # its.
'subspaceMin',      0 , ... % Use subspace minimization
'iscomplex'  ,    NaN , ... % Flag set to indicate complex problem
'maxMatvec'  ,    Inf , ... % Maximum matrix-vector multiplies allowed
'weights'    ,      1 , ... % Weights W in ||Wx||_1
'project'    , @NormL1_project , ...
'primal_norm', @NormL1_primal  , ...
'dual_norm'  , @NormL1_dual      ...
   );
options = spgSetParms(defaultopts, options);

fid           = options.fid;
logLevel      = options.verbosity;
maxIts        = options.iterations;
nPrevVals     = options.nPrevVals;
bpTol         = options.bpTol;
lsTol         = options.lsTol;
optTol        = options.optTol;
decTol        = options.decTol;
stepMin       = options.stepMin;
stepMax       = options.stepMax;
activeSetIt   = options.activeSetIt;
subspaceMin   = options.subspaceMin;
maxMatvec     = max(3,options.maxMatvec);
weights       = options.weights;

maxLineErrors = 10;     % Maximum number of line-search failures.
pivTol        = 1e-12;  % Threshold for significant Newton step.

%----------------------------------------------------------------------
% Initialize local variables.
%----------------------------------------------------------------------
iter          = 0;  itnTotLSQR = 0; % Total SPGL1 and LSQR iterations.
nProdA        = 0;  nProdAt    = 0;
lastFv        = -inf(nPrevVals,1);  % Last m function values.
nLineTot      = 0;                  % Total no. of linesearch steps.
printTau      = false;
nNewton       = 0;
bNorm         = norm(b,2);
stat          = false;
timeProject   = 0;
timeMatProd   = 0;
nnzIter       = 0;                  % No. of its with fixed pattern.
nnzIdx        = [];                 % Active-set indicator.
subspace      = false;              % Flag if did subspace min in current itn.
stepG         = 1;                  % Step length for projected gradient.
testUpdateTau = 0;                  % Previous step did not update tau

% Determine initial x, vector length n, and see if problem is complex
explicit = ~(isa(A,'function_handle'));
if isempty(x)
   if isnumeric(A)
      n = size(A,2);
      realx = isreal(A) && isreal(b);
   else
      x = Aprod(b,2);
      n = length(x);
      realx = isreal(x) && isreal(b);
   end
   x = zeros(n,1);
else
   n     = length(x);
   realx = isreal(x) && isreal(b);
end
if isnumeric(A), realx = realx && isreal(A); end;

% Override options when options.iscomplex flag is set
if (~isnan(options.iscomplex)), realx = (options.iscomplex == 0); end

% Check if all weights (if any) are strictly positive. In previous
% versions we also checked if the number of weights was equal to
% n. In the case of multiple measurement vectors, this no longer
% needs to apply, so the check was removed.
if ~isempty(weights)
  if any(~isfinite(weights))
     error('Entries in options.weights must be finite');
  end
  if any(weights <= 0)
     error('Entries in options.weights must be strictly positive');
  end
else
  weights = 1;
end

% Quick exit if sigma >= ||b||.  Set tau = 0 to short-circuit the loop.
if bNorm <= sigma
   printf('W: sigma >= ||b||.  Exact solution is x = 0.\n');
   tau = 0;  singleTau = true;
end 
  
% Don't do subspace minimization if x is complex.
if ~realx && subspaceMin
   printf('W: Subspace minimization disabled when variables are complex.\n');
   subspaceMin = false;
end

% Pre-allocate iteration info vectors
xNorm1 = zeros(min(maxIts,10000),1);
rNorm2 = zeros(min(maxIts,10000),1);
lambda = zeros(min(maxIts,10000),1);

% Exit conditions (constants).
EXIT_ROOT_FOUND    = 1;
EXIT_BPSOL_FOUND   = 2;
EXIT_LEAST_SQUARES = 3;
EXIT_OPTIMAL       = 4;
EXIT_ITERATIONS    = 5;
EXIT_LINE_ERROR    = 6;
EXIT_SUBOPTIMAL_BP = 7;
EXIT_MATVEC_LIMIT  = 8;
EXIT_ACTIVE_SET    = 9; % [sic]

%----------------------------------------------------------------------
% Log header.
%----------------------------------------------------------------------
printf('\n');
printf(' %s\n',repmat('=',1,80));
printf(' SPGL1  v.%s (%s)\n', REVISION, DATE);
printf(' %s\n',repmat('=',1,80));
printf(' %-22s: %8i %4s'   ,'No. rows'          ,m       ,'');
printf(' %-22s: %8i\n'     ,'No. columns'       ,n          );
printf(' %-22s: %8.2e %4s' ,'Initial tau'       ,tau     ,'');
printf(' %-22s: %8.2e\n'   ,'Two-norm of b'     ,bNorm      );
printf(' %-22s: %8.2e %4s' ,'Optimality tol'    ,optTol  ,'');
if singleTau
   printf(' %-22s: %8.2e\n'  ,'Target one-norm of x'  ,tau       );
else
   printf(' %-22s: %8.2e\n','Target objective'  ,sigma      );
end
printf(' %-22s: %8.2e %4s' ,'Basis pursuit tol' ,bpTol   ,'');
printf(' %-22s: %8i\n'     ,'Maximum iterations',maxIts     );
printf('\n');
if singleTau
   logB = ' %5i  %13.7e  %13.7e  %9.2e  %6.1f  %6i  %6i';
   logH = ' %5s  %13s  %13s  %9s  %6s  %6s  %6s\n';
   printf(logH,'Iter','Objective','Relative Gap','gNorm','stepG','nnzX','nnzG');
else
   logB = ' %5i  %13.7e  %13.7e  %9.2e  %9.3e  %6.1f  %6i  %6i';
   logH = ' %5s  %13s  %13s  %9s  %9s  %6s  %6s  %6s  %13s\n';
   printf(logH,'Iter','Objective','Relative Gap','Rel Error',...
          'gNorm','stepG','nnzX','nnzG','tau');
end    
    
% Project the starting point and evaluate function and gradient.
x         = project(x,tau);
r         = b - Aprod(x,1);  % r = b - Ax
g         =   - Aprod(r,2);  % g = -A'r
f         = r'*r / 2; 

% Required for nonmonotone strategy.
lastFv(1) = f;
fBest     = f;
xBest     = x;
fOld      = f;

% Compute projected gradient direction and initial steplength.
dx     = project(x - g, tau) - x;
dxNorm = norm(dx,inf);
if dxNorm < (1 / stepMax)
   gStep = stepMax;
else
   gStep = min( stepMax, max(stepMin, 1/dxNorm) );
end
