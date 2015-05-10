
%----------------------------------------------------------------------
% MAIN LOOP.
%----------------------------------------------------------------------
while 1
    
    %------------------------------------------------------------------
    % Test exit conditions.
    %------------------------------------------------------------------

    % Compute quantities needed for log and exit conditions.
    gNorm   = options.dual_norm(-g,weights);
    rNorm   = norm(r, 2);
    gap     = r'*(r-b) + tau*gNorm;
    rGap    = abs(gap) / max(1,f);
    aError1 = rNorm - sigma;
    aError2 = f - sigma^2 / 2;
    rError1 = abs(aError1) / max(1,rNorm);
    rError2 = abs(aError2) / max(1,f);

    % Count number of consecutive iterations with identical support.
    nnzOld = nnzIdx;
    [nnzX,nnzG,nnzIdx,nnzDiff] = activeVars(x,g,nnzIdx,options);

    if nnzDiff
       nnzIter = 0;
    else
       nnzIter = nnzIter + 1;
       if nnzIter >= activeSetIt, stat=EXIT_ACTIVE_SET; end
    end
    
    % Single tau: Check if we're optimal.
    % The 2nd condition is there to guard against large tau.
    if singleTau
       if rGap <= optTol || rNorm < optTol*bNorm
          stat  = EXIT_OPTIMAL;
       end
 
    % Multiple tau: Check if found root and/or if tau needs updating.
    else
       
       % Test if a least-squares solution has been found
       if gNorm <= lsTol * rNorm
          stat = EXIT_LEAST_SQUARES;
       end
       
       if rGap <= max(optTol, rError2) || rError1 <= optTol
          % The problem is nearly optimal for the current tau.
          % Check optimality of the current root.
          test1 = rNorm       <=   bpTol * bNorm;
       %  test2 = gNorm       <=   bpTol * rNorm;
          test3 = rError1     <=  optTol;
          test4 = rNorm       <=  sigma;
          
          if test4, stat=EXIT_SUBOPTIMAL_BP;end  % Found suboptimal BP sol.
          if test3, stat=EXIT_ROOT_FOUND;   end  % Found approx root.
          if test1, stat=EXIT_BPSOL_FOUND;  end  % Resid minim'zd -> BP sol.
        % 30 Jun 09: Large tau could mean large rGap even near LS sol.
        %            Move LS check out of this if statement.
        % if test2, stat=EXIT_LEAST_SQUARES; end % Gradient zero -> BP sol.
       end
       
       testRelChange1 = (abs(f - fOld) <= decTol * f);
       testRelChange2 = (abs(f - fOld) <= 1e-1 * f * (abs(rNorm - sigma)));
       testUpdateTau  = ((testRelChange1 && rNorm >  2 * sigma) || ...
                         (testRelChange2 && rNorm <= 2 * sigma)) && ...
                         ~stat && ~testUpdateTau;
       
       if testUpdateTau
          % Update tau.
          tauOld   = tau;
          tau      = max(0,tau + (rNorm * aError1) / gNorm);
          nNewton  = nNewton + 1;
          printTau = abs(tauOld - tau) >= 1e-6 * tau; % For log only.
          if tau < tauOld
             % The one-norm ball has decreased.  Need to make sure that the
             % next iterate if feasible, which we do by projecting it.
             x = project(x,tau);
             
             % Update the residual, gradient, and function value.
             r = b - Aprod(x,1);  % r = b - Ax
             g =   - Aprod(r,2);  % g = -A'r
             f = r'*r / 2;
             
             % Reset the function value history.
             lastFv = -inf(nPrevVals,1);  % Last m function values.
             lastFv(1) = f;
          end
       end
    end

    % Too many its and not converged.
    if ~stat  &&  iter >= maxIts
        stat = EXIT_ITERATIONS;
    end

    %------------------------------------------------------------------
    % Print log, update history and act on exit conditions.
    %------------------------------------------------------------------
    if logLevel >= 2 || singleTau || printTau || iter == 0 || stat
       tauFlag = '              '; subFlag = '';
       if printTau, tauFlag = sprintf(' %13.7e',tau);   end
       if subspace, subFlag = sprintf(' S %2i',itnLSQR); end
       if singleTau
          printf(logB,iter,rNorm,rGap,gNorm,log10(stepG),nnzX,nnzG);
          if subspace
             printf('  %s',subFlag);
          end
       else
          printf(logB,iter,rNorm,rGap,rError1,gNorm,log10(stepG),nnzX,nnzG);
          if printTau || subspace
             printf(' %s',[tauFlag subFlag]);
          end
       end
       printf('\n');
    end
    printTau = false;
    subspace = false;
    
    % Update history info
    xNorm1(iter+1) = options.primal_norm(x,weights);
    rNorm2(iter+1) = rNorm;
    lambda(iter+1) = gNorm;
    
    if stat, break; end % Act on exit conditions.
        
    %==================================================================
    % Iterations begin here.
    %==================================================================
    iter = iter + 1;
    xOld = x;  fOld = f;  gOld = g;  rOld = r;

    try
       %---------------------------------------------------------------
       % Projected gradient step and linesearch.
       %---------------------------------------------------------------
       [f,x,r,nLine,stepG,lnErr] = ...
           spgLineCurvy(x,gStep*g,max(lastFv),@Aprod,b,@project,tau);
       nLineTot = nLineTot + nLine;
       if lnErr
          % Projected backtrack failed. Retry with feasible dir'n linesearch.
          x    = xOld;
          f    = fOld;
          r    = rOld;
          dx   = project(x - gStep*g, tau) - x;
          gtd  = g'*dx;
          [f,x,r,nLine,lnErr] = spgLine(f,x,dx,gtd,max(lastFv),@Aprod,b);
          nLineTot = nLineTot + nLine;
       end
       if lnErr
          % Failed again. Revert to previous iterates and damp max BB step.
          x = xOld;
          f = fOld;
          if maxLineErrors <= 0
             stat = EXIT_LINE_ERROR;
          else
             stepMax = stepMax / 10;
             printf(['W: Linesearch failed with error %i. '...
                     'Damping max BB scaling to %6.1e.\n'],lnErr,stepMax);
             maxLineErrors = maxLineErrors - 1;
          end
       end

       %---------------------------------------------------------------
       % Subspace minimization (only if active-set change is small).
       %---------------------------------------------------------------
       doSubspaceMin = false;
       if subspaceMin
          g = - Aprod(r,2);
          [nnzX,nnzG,nnzIdx,nnzDiff] = activeVars(x,g,nnzOld,options);
          if ~nnzDiff
              if nnzX == nnzG, itnMaxLSQR = 20;
              else             itnMaxLSQR = 5;
              end
              nnzIdx = abs(x) >= optTol; 
              doSubspaceMin = true;
          end
       end

       if doSubspaceMin
     
          % LSQR parameters
          damp       = 1e-5;
          aTol       = 1e-1;
          bTol       = 1e-1;
          conLim     = 1e12;
          showLSQR   = 0;
       
          ebar   = sign(x(nnzIdx));
          nebar  = length(ebar);
          Sprod  = @(y,mode)LSQRprod(@Aprod,nnzIdx,ebar,n,y,mode);
       
          [dxbar, istop, itnLSQR] = ...
             lsqr(m,nebar,Sprod,r,damp,aTol,bTol,conLim,itnMaxLSQR,showLSQR);
              
          itnTotLSQR = itnTotLSQR + itnLSQR;
       
          if istop ~= 4  % LSQR iterations successful. Take the subspace step.
             % Push dx back into full space:  dx = Z dx.
             dx = zeros(n,1);
             dx(nnzIdx) = dxbar - (1/nebar)*(ebar'*dxbar)*ebar;

             % Find largest step to a change in sign.
             block1 = nnzIdx  &  x < 0  &  dx > +pivTol;
             block2 = nnzIdx  &  x > 0  &  dx < -pivTol;
             alpha1 = Inf; alpha2 = Inf;
             if any(block1), alpha1 = min(-x(block1) ./ dx(block1)); end
             if any(block2), alpha2 = min(-x(block2) ./ dx(block2)); end
             alpha = min([1  alpha1  alpha2]);
             ensure(alpha >= 0);
             ensure(ebar'*dx(nnzIdx) <= optTol);          
          
             % Update variables.
             x    = x + alpha*dx;
             r    = b - Aprod(x,1);
             f    = r'*r / 2;
             subspace = true;
          end
       end
       
       ensure(options.primal_norm(x,weights) <= tau+optTol);

       %---------------------------------------------------------------
       % Update gradient and compute new Barzilai-Borwein scaling.
       %---------------------------------------------------------------
       if (~lnErr)
          g    = - Aprod(r,2);
          s    = x - xOld;
          y    = g - gOld;
          sts  = s'*s;
          sty  = s'*y;
          if   sty <= 0,  gStep = stepMax;
          else            gStep = min( stepMax, max(stepMin, sts/sty) );
          end
       else
          gStep = min( stepMax, gStep );
       end
       
    catch err % Detect matrix-vector multiply limit error
       if strcmp(err.identifier,'SPGL1:MaximumMatvec')
         stat = EXIT_MATVEC_LIMIT;
         iter = iter - 1;
         x = xOld;  f = fOld;  g = gOld;  r = rOld;
         break;
       else
         rethrow(err);
       end
    end

    %------------------------------------------------------------------
    % Update function history.
    %------------------------------------------------------------------
    if singleTau || f > sigma^2 / 2 % Don't update if superoptimal.
       lastFv(mod(iter,nPrevVals)+1) = f;
       if fBest > f
          fBest = f;
          xBest = x;
       end
    end

end % while 1