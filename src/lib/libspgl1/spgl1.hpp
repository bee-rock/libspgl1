#include <algorithm>

namespace libspgl1 {


template<typename MatrixType, typename VectorType>
VectorType spgl1(const MatrixType& A, const MatrixType& At, const VectorType& b, const VectorType& x0){
	libspgl1::Parameters parameters = libspgl1::Parameters(b);
	VectorType sign_x = libspgl1::vector::sign<VectorType>(x0);
	VectorType x = libspgl1::projectI<VectorType>(x0, parameters.tau);
	VectorType r = libspgl1::initialization::compute_r(A, b, x);
	double f = libspgl1::initialization::compute_f(r);
	VectorType g = libspgl1::initialization::compute_g(At, r);
	VectorType dx = libspgl1::projectI<VectorType>(x-g, parameters.tau)-x;
	double dxNorm = libspgl1::math::max<double>(libspgl1::vector::abs<VectorType>(dx));
	double gStep = libspgl1::initialization::compute_gstep(dxNorm,
														   parameters.stepmin,
														   parameters.stepmax);
	std::cout << "f: " << f << std::endl;
	std::cout << "gstep: " << gStep << std::endl;
	double lastFv    = f;
	double fBest     = f;
	VectorType xBest = x;
	double fOld      = f;
	return dx;
}

} // libspgl1
