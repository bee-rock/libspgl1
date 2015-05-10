
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS.  These share some vars with workspace above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function z = Aprod(x,mode)
   if (nProdA + nProdAt >= maxMatvec)
     error('SPGL1:MaximumMatvec','');
   end
     
   tStart = toc;
   if mode == 1
      nProdA = nProdA + 1;
      if   explicit, z = A*x;
      else           z = A(x,1);
      end
   elseif mode == 2
      nProdAt = nProdAt + 1;
      if   explicit, z = A'*x;
      else           z = A(x,2);
      end
   else
      error('Wrong mode!');
   end
   timeMatProd = timeMatProd + (toc - tStart);
end % function Aprod

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printf(varargin)
  if logLevel > 0
     fprintf(fid,varargin{:});
  end
end % function printf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = project(x, tau)
   tStart      = toc;

   x = options.project(x,weights,tau);
   
   timeProject = timeProject + (toc - tStart);
end % function project

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of nested functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % function spg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nnzX,nnzG,nnzIdx,nnzDiff] = activeVars(x,g,nnzIdx,options)
% Find the current active set.
% nnzX    is the number of nonzero x.
% nnzG    is the number of elements in nnzIdx.
% nnzIdx  is a vector of primal/dual indicators.
% nnzDiff is the no. of elements that changed in the support.
  xTol    = min(.1,10*options.optTol);
  gTol    = min(.1,10*options.optTol);
  gNorm   = options.dual_norm(g,options.weights);
  nnzOld  = nnzIdx;

  % Reduced costs for postive & negative parts of x.
  z1 = gNorm + g;
  z2 = gNorm - g;

  % Primal/dual based indicators.
  xPos    = x >  xTol  &  z1 < gTol; %g < gTol;%
  xNeg    = x < -xTol  &  z2 < gTol; %g > gTol;%
  nnzIdx  = xPos | xNeg;

  % Count is based on simple primal indicator.
  nnzX    = sum(abs(x) >= xTol);
  nnzG    = sum(nnzIdx);
  
  if isempty(nnzOld)
     nnzDiff = inf;
  else
     nnzDiff = sum(nnzIdx ~= nnzOld);
  end
  
end % function spgActiveVars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = LSQRprod(Aprod,nnzIdx,ebar,n,dx,mode)
% Matrix multiplication for subspace minimization.
% Only called by LSQR.
  nbar = length(ebar);
   if mode == 1
      y = zeros(n,1);
      y(nnzIdx) = dx - (1/nbar)*(ebar'*dx)*ebar; % y(nnzIdx) = Z*dx
      z = Aprod(y,1);                            % z = S Z dx
   else
      y = Aprod(dx,2);
      z = y(nnzIdx) - (1/nbar)*(ebar'*y(nnzIdx))*ebar;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fNew,xNew,rNew,iter,err] = spgLine(f,x,d,gtd,fMax,Aprod,b)
% Nonmonotone linesearch.

EXIT_CONVERGED  = 0;
EXIT_ITERATIONS = 1;
maxIts = 10;
step   = 1;
iter   = 0;
gamma  = 1e-4;
gtd    = -abs(gtd); % 03 Aug 07: If gtd is complex,
                    % then should be looking at -abs(gtd).
while 1

    % Evaluate trial point and function value.
    xNew = x + step*d;
    rNew = b - Aprod(xNew,1);
    fNew = rNew'*rNew / 2;

    % Check exit conditions.
    if fNew < fMax + gamma*step*gtd  % Sufficient descent condition.
       err = EXIT_CONVERGED;
       break
    elseif  iter >= maxIts           % Too many linesearch iterations.
       err = EXIT_ITERATIONS;
       break
    end
    
    % New linesearch iteration.
    iter = iter + 1;
    
    % Safeguarded quadratic interpolation.
    if step <= 0.1
       step  = step / 2;
    else
       tmp = (-gtd*step^2) / (2*(fNew-f-step*gtd));
       if tmp < 0.1 || tmp > 0.9*step || isnan(tmp)
          tmp = step / 2;
       end
       step = tmp;
    end
    
end % while 1

end % function spgLine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fNew,xNew,rNew,iter,step,err] = ...
    spgLineCurvy(x,g,fMax,Aprod,b,project,tau)
% Projected backtracking linesearch.
% On entry,
% g  is the (possibly scaled) steepest descent direction.

EXIT_CONVERGED  = 0;
EXIT_ITERATIONS = 1;
EXIT_NODESCENT  = 2;
gamma  = 1e-4;
maxIts = 10;
step   =  1;
sNorm  =  0;
scale  =  1;      % Safeguard scaling.  (See below.)
nSafe  =  0;      % No. of safeguarding steps.
iter   =  0;
debug  =  false;  % Set to true to enable log.
n      =  length(x);

if debug
   fprintf(' %5s  %13s  %13s  %13s  %8s\n',...
           'LSits','fNew','step','gts','scale');  
end
   
while 1

    % Evaluate trial point and function value.
    xNew     = project(x - step*scale*g, tau);
    rNew     = b - Aprod(xNew,1);
    fNew     = rNew'*rNew / 2;
    s        = xNew - x;
    gts      = scale * real(g' * s);
%   gts      = scale * (g' * s);
    if gts >= 0
       err = EXIT_NODESCENT;
       break
    end

    if debug
       fprintf(' LS %2i  %13.7e  %13.7e  %13.6e  %8.1e\n',...
               iter,fNew,step,gts,scale);
    end
    
    % 03 Aug 07: If gts is complex, then should be looking at -abs(gts).
    % 13 Jul 11: It's enough to use real part of g's (see above).
    if fNew < fMax + gamma*step*gts
%   if fNew < fMax - gamma*step*abs(gts)  % Sufficient descent condition.
       err = EXIT_CONVERGED;
       break
    elseif iter >= maxIts                 % Too many linesearch iterations.
       err = EXIT_ITERATIONS;
       break
    end
    
    % New linesearch iteration.
    iter = iter + 1;
    step = step / 2;

    % Safeguard: If stepMax is huge, then even damped search
    % directions can give exactly the same point after projection.  If
    % we observe this in adjacent iterations, we drastically damp the
    % next search direction.
    % 31 May 07: Damp consecutive safeguarding steps.
    sNormOld  = sNorm;
    sNorm     = norm(s) / sqrt(n);
    %   if sNorm >= sNormOld
    if abs(sNorm - sNormOld) <= 1e-6 * sNorm
       gNorm = norm(g) / sqrt(n);
       scale = sNorm / gNorm / (2^nSafe);
       nSafe = nSafe + 1;
    end
    
end % while 1

end % function spgLineCurvy


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = NormL1_project(x,weights,tau)

if isreal(x)
   x = oneProjector(x,weights,tau);
else
   xa  = abs(x);
   idx = xa < eps;
   xc  = oneProjector(xa,weights,tau);
   xc  = xc ./ xa; xc(idx) = 0;
   x   = x .* xc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = NormL1_primal(x,weights)

p = norm(x.*weights,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = NormL1_project(x,weights,tau)

if isreal(x)
   x = oneProjector(x,weights,tau);
else
   xa  = abs(x);
   idx = xa < eps;
   xc  = oneProjector(xa,weights,tau);
   xc  = xc ./ xa; xc(idx) = 0;
   x   = x .* xc;
end
