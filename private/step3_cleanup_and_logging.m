
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Restore best solution (only if solving single problem).
if singleTau && f > fBest
   rNorm = sqrt(2*fBest);
   printf('\n Restoring best iterate to objective %13.7e\n',rNorm);
   x = xBest;
   r = b - Aprod(x,1);
   g =   - Aprod(r,2);
   gNorm = options.dual_norm(g,weights);
   rNorm = norm(r,  2);
end

% Final cleanup before exit.
info.tau         = tau;
info.rNorm       = rNorm;
info.rGap        = rGap;
info.gNorm       = gNorm;
info.rGap        = rGap;
info.stat        = stat;
info.iter        = iter;
info.nProdA      = nProdA;
info.nProdAt     = nProdAt;
info.nNewton     = nNewton;
info.timeProject = timeProject;
info.timeMatProd = timeMatProd;
info.itnLSQR     = itnTotLSQR;
info.options     = options;
info.timeTotal   = toc;

info.xNorm1      = xNorm1(1:iter);
info.rNorm2      = rNorm2(1:iter);
info.lambda      = lambda(1:iter);

% Print final output.
switch (stat)
   case EXIT_OPTIMAL
      printf('\n EXIT -- Optimal solution found\n')
   case EXIT_ITERATIONS
      printf('\n ERROR EXIT -- Too many iterations\n');
   case EXIT_ROOT_FOUND
      printf('\n EXIT -- Found a root\n');
   case {EXIT_BPSOL_FOUND}
      printf('\n EXIT -- Found a BP solution\n');
   case {EXIT_LEAST_SQUARES}
      printf('\n EXIT -- Found a least-squares solution\n');
   case EXIT_LINE_ERROR
      printf('\n ERROR EXIT -- Linesearch error (%i)\n',lnErr);
   case EXIT_SUBOPTIMAL_BP
      printf('\n EXIT -- Found a suboptimal BP solution\n');
   case EXIT_MATVEC_LIMIT
      printf('\n EXIT -- Maximum matrix-vector operations reached\n');
   case EXIT_ACTIVE_SET
      printf('\n EXIT -- Found a possible active set\n');
   otherwise
      error('Unknown termination condition\n');
end
printf('\n');
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Products with A',nProdA,'','Total time   (secs)',info.timeTotal);
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Products with A''',nProdAt,'','Project time (secs)',timeProject);
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Newton iterations',nNewton,'','Mat-vec time (secs)',timeMatProd);
printf(' %-20s:  %6i %6s %-20s:  %6i\n', ...
   'Line search its',nLineTot,'','Subspace iterations',itnTotLSQR);
printf('\n');


