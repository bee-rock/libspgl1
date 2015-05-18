#pragma once
#include "vector_adapter.hpp"

namespace libspgl1 {
namespace math {

template<typename ElementType>
ElementType pow(const ElementType& x, const ElementType& p) {
    return std::pow(x, p);
}

template<typename ElementType, typename VectorType>
double norm(const VectorType& a, const double& p){
	double result{0};
	const size_t n_elems = libspgl1::vector::n_elem(a);
	for (size_t i=0; i<n_elems; ++i)
		result += libspgl1::math::pow<ElementType>(
				std::abs(libspgl1::vector::get_element<ElementType>(a, i)),p);
	return static_cast<ElementType>(libspgl1::math::pow<ElementType>(result, 1/p));
}

template<typename ElementType, typename VectorType>
ElementType NormL1_primal(const VectorType& a, const VectorType& b)
{
	size_t n_elems = libspgl1::vector::n_elem(a);
	ElementType result{0};
	for (size_t i=0; i<n_elems; ++i)
		result += std::abs(libspgl1::vector::get_element<ElementType>(a, i)) * libspgl1::vector::get_element<ElementType>(b, i);
	return result;
}

} // math
} // libspgl1
