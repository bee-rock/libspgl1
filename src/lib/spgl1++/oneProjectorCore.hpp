#pragma once
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <sys/types.h>
#include <algorithm>

namespace libspgl1 {

template<typename VectorType>
double calculate_threshhold_value(const VectorType &c_bar, double tau){
	double csb = 0.0 - tau;
	double soft_threshold_value = 0;
	double alpha = 0.0;
	for (size_t j = 0; j < libspgl1::vector::n_elem(c_bar); j++)
	{
		/* Get current maximum heap element         */
		double b = libspgl1::vector::get_element<double>(c_bar, libspgl1::vector::n_elem(c_bar)-1-j);
		csb += b;          /* Update the cumulative sum of b           */

		/* Compute the required step to satisfy the tau constraint */
		alpha  = csb / (static_cast<double>(j)+1.0);

		/* We are done as soon as the constraint can be satisfied    *
		 * /* without exceeding the current minimum value in `vector' b */
		if (alpha >= b){
			return soft_threshold_value;
		}
		soft_threshold_value = alpha;
	}
	return soft_threshold_value;
}

template<typename VectorType>
VectorType projectI(const VectorType& c, const double tau)
{
	double csb  = 0.0;        /* Cumulative sum of b */
	double alpha = 0.0;
	VectorType sign_c = libspgl1::vector::sign<VectorType>(c);
	VectorType c_bar = libspgl1::vector::abs(c);
	size_t n = libspgl1::vector::n_elem(c);

   /* Check if tau is essentially zero.  Exit with x = 0. */
	if (tau < DBL_EPSILON) {
		for (size_t i = 0; i < n; i++)
		{
			libspgl1::vector::set_element<double>(c_bar, i, 0.0);
		}
		return c_bar;
	}

   /* Check if ||b||_1 <= lambda.  Exit with x = b. */
	for (size_t i = 0; i < n; i++){
		csb += libspgl1::vector::get_element<double>(c_bar, i);
	}
	if (csb <= tau){
		return c;
	}

	std::make_heap(c_bar.begin(), c_bar.end());
	std::sort_heap(c_bar.begin(), c_bar.end());

	double soft_threshold_value = calculate_threshhold_value<VectorType>(c_bar, tau);

	/* Set the solution by applying soft-thresholding with `soft' */
	for (size_t i = 0; i < n; i++)
	{
		double b = std::abs(libspgl1::vector::get_element<double>(c, i));
		if (b <= soft_threshold_value){
			libspgl1::vector::set_element<double>(c_bar, i, 0.0);
		}
		else{
			libspgl1::vector::set_element<double>(c_bar, i, b - soft_threshold_value);
		}
	}
	return libspgl1::vector::elementwise_multiplication<VectorType>(c_bar, sign_c);
}

} // libspgl1
