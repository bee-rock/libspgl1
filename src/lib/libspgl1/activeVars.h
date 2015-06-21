//function [,,,] = activeVars(x,g,nnzIdx,options)

namespace libspgl1 {

	class ActiveVars {
		size_t nnzOld, nnzX, nnzG, nnzIdx, nnzDiff;

		ActiveVars() :
			nnzOld{0}, nnzX{0}, nnzG{0}, nnzIdx{0}, nnzDiff{0}
		{}

		size_t compute_nnzX();
		size_t compute_nnzG();
		size_t compute_nnzIdx();
		size_t compute_nnzDiff();
		double xTol  = std::min(.1, 10.0*options.optTol);
		double  gTol = std::min(.1, 10.0*options.optTol);
		double  gNorm   = options.dual_norm(g, options.weights);
		nnzOld = nnzIdx;
		//  % Reduced costs for postive & negative parts of x.
		double  z1 = gNorm + g;
		double  z2 = gNorm - g;
		//
		//  % Primal/dual based indicators.
		//  xPos    = x >  xTol  &  z1 < gTol; %g < gTol;%
		//  xNeg    = x < -xTol  &  z2 < gTol; %g > gTol;%
		//  nnzIdx  = xPos | xNeg;
		//
		//  % Count is based on simple primal indicator.
		//  nnzX    = sum(abs(x) >= xTol);
		//  nnzG    = sum(nnzIdx);
		//
		//  if isempty(nnzOld)
		//     nnzDiff = inf;
		//  else
		//     nnzDiff = sum(nnzIdx ~= nnzOld);
		//  end
	};

} // libspgl1
